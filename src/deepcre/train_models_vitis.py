from __future__ import annotations
import argparse
import json
import os
import tensorflow as tf                                                                                          #type:ignore
from typing import Any, Dict, List, Optional, Tuple
import pandas as pd
from tensorflow.keras.layers import Dropout, Dense, Input, Conv1D, Activation, MaxPool1D, Flatten               #type:ignore
from tensorflow.keras.optimizers import Adam                                                                    #type:ignore
from tensorflow.keras import Model                                                                              #type:ignore
from tensorflow.keras.models import load_model                                                                  #type:ignore
from tensorflow.keras.callbacks import ModelCheckpoint, EarlyStopping, ReduceLROnPlateau                        #type:ignore
from tensorflow.keras.metrics import AUC     
import gc
from tensorflow.keras import backend as K                                                                   #type:ignore
import pickle
import numpy as np
from pyfaidx import Fasta
from sklearn.utils import shuffle
import re
import sys

from parsing import ParsedInputs, RunInfo, ModelCase
from utils import get_filename_from_path, get_time_stamp, one_hot_encode, make_absolute_path, result_summary, combine_annotations, combine_fasta, combine_tpms, load_input_files

# Global tracking for consistently excluded genes across all models
GLOBAL_EXCLUDED_GENES = {
    'consistently_excluded': set(),  # truely excluded genes (beyond expected validation exclusions)
    'model_count': 0,  # Number of models trained
    'timestamp': None
}

def configure_gpu_memory():
    """Configure GPU to avoid memory accumulation"""
    gpus = tf.config.experimental.list_physical_devices('GPU')
    if gpus:
        try:
            # Enable memory growth
            for gpu in gpus:
                tf.config.experimental.set_memory_growth(gpu, True)
        except RuntimeError as e:
            print(f"GPU configuration error: {e}")

def cleanup_memory():
    """Memory cleanup between chromosome training"""
    K.clear_session()  # Clear TensorFlow session
    gc.collect()       # Force garbage collection
    
    try:
        tf.keras.backend.clear_session()
        if tf.config.experimental.list_physical_devices('GPU'):
            tf.config.experimental.reset_memory_stats('GPU:0')
    except:
        pass

def find_completed_chromosomes(output_name: str, model_case: ModelCase) -> List[str]:
    """Find which chromosomes already have trained models"""
    models_path = make_absolute_path("../../out/saved_models", start_file=__file__)
    completed = []
    
    if not os.path.exists(models_path):
        return completed
    
    for file in os.listdir(models_path):
        if output_name in file and str(model_case) in file:
            # Extract chromosome from filename
            match = re.search(rf'{output_name}_(.+?)_{model_case}_', file)
            if match:
                completed.append(match.group(1))
    
    return completed

def filter_remaining_chromosomes(chromosomes: List[str], output_name: str, model_case: ModelCase) -> List[str]:
    """Filter out already completed chromosomes"""
    completed = find_completed_chromosomes(output_name, model_case)
    remaining = [chr for chr in chromosomes if chr not in completed]
    
    print(f"Found {len(completed)} completed chromosomes: {completed}")
    print(f"Remaining chromosomes to train: {remaining}")
    
    return remaining

class TerminationError(Exception):
    """Error that warrants termination of the whole program"""
    pass

def find_newest_model_path(output_name: str, model_case: ModelCase, val_chromosome: str = "", model_path: str = "") -> Dict[str, str]:
    """finds path to newest model fitting the given parameters"""
    if model_path == "":
        path_to_models = make_absolute_path("../../out/saved_models", start_file=__file__)
    else:
        path_to_models = make_absolute_path(model_path, start_file=__file__)
    
    if val_chromosome == "":
        regex_string = f"^{output_name}_(.+)_{model_case}_train_models_vitis_\d+_\d+\.h5$"
    else:
        regex_string = f"^{output_name}_{val_chromosome}_{model_case}_train_models_vitis_\d+_\d+\.h5$"
        
    regex = re.compile(regex_string)
    candidate_models = [model for model in os.listdir(path_to_models)]
    fitting_models = {}
    
    for candidate in candidate_models:
        match = regex.match(candidate)
        if match:
            if val_chromosome:
                chromosome = val_chromosome
            else:
                chromosome = match.group(1)
            if chromosome in fitting_models:
                fitting_models[chromosome].append(candidate)
            else:
                fitting_models[chromosome] = [candidate]

    if not fitting_models:
        raise ValueError(f"no trained models fitting the given parameters (output_name: '{output_name}', val_chromosome: '{val_chromosome}', model_case: '{model_case}') were found!")
    
    for chromosome, models in fitting_models.items():
        models.sort()
        fitting_models[chromosome] = os.path.join(path_to_models, models[-1])
    return fitting_models

def extract_gene(genome: Fasta, extragenic: int, intragenic: int, ignore_small_genes: bool, expected_final_size: int,
                 chrom: str, start: int, end: int, strand: str, ssc_training: bool = False, val_chromosome: Optional[str] = None) -> np.ndarray:
    """extracts the gene flanking region for a single gene and converts it to a numpy encoded one-hot encoding"""
    gene_size = end - start
    extractable_intragenic = intragenic if gene_size // 2 > intragenic else gene_size // 2
    prom_start, prom_end = start - extragenic, start + extractable_intragenic
    term_start, term_end = end - extractable_intragenic, end + extragenic

    promoter = one_hot_encode(genome[chrom][prom_start:prom_end])
    terminator = one_hot_encode(genome[chrom][term_start:term_end])
    extracted_size = promoter.shape[0] + terminator.shape[0]
    central_pad_size = expected_final_size - extracted_size

    if ssc_training and chrom != val_chromosome:
            np.random.shuffle(promoter)
            np.random.shuffle(terminator)

    pad_size = 20 if ignore_small_genes else central_pad_size

    if strand == '+':
        seq = np.concatenate([
            promoter,
            np.zeros(shape=(pad_size, 4)),
            terminator
        ])
    else:
        seq = np.concatenate([
            terminator[::-1, ::-1],
            np.zeros(shape=(pad_size, 4)),
            promoter[::-1, ::-1]
        ])
        
    return seq

def append_sequence_training(include_as_validation_gene: bool, include_as_training_gene: bool, expected_final_size: int,
                             train_seqs: List, val_seqs: List, train_targets: List, val_targets: List, tpms: pd.DataFrame,
                             gene_id: str, seq: np.ndarray) -> Tuple[int, int]:
    """appends a sequence to the training or validation set if it has the expected size"""
    added_val, added_training = 0, 0
    if seq.shape[0] == expected_final_size:
        if include_as_validation_gene:
            val_seqs.append(seq)
            val_targets.append(tpms.loc[gene_id, 'target'])
            added_val = 1
        elif include_as_training_gene:
            train_seqs.append(seq)
            train_targets.append(tpms.loc[gene_id, 'target'])
            added_training = 1
    return added_val, added_training

def calculate_conditions(val_chromosome: str, model_case: ModelCase, train_val_split: bool, val_specie: str,
                         validation_genes: List, current_val_size: int, current_train_size: int, target_val_size: int,
                         target_train_size: int, specie: str, chrom: str, gene_id: str) -> Tuple[bool, bool]:
    """calculates whether a gene should be included in the validation or training set"""
    if model_case == ModelCase.MSR:
        include_in_validation_set = specie == val_specie
        include_in_training_set = not include_in_validation_set
    elif not train_val_split:
        # ORIGINAL CHROMOSOME-BASED LOGIC - NO SIZE LIMITS
        include_in_validation_set = chrom == val_chromosome and gene_id in validation_genes
        include_in_training_set = chrom != val_chromosome
    else:
        # Random split with size limits (kept for when train_val_split=True)
        include_in_validation_set = current_val_size < target_val_size and gene_id in validation_genes
        include_in_training_set = not include_in_validation_set and current_train_size < target_train_size
    
    return include_in_validation_set, include_in_training_set

def set_up_validation_genes(pickled_genes_path: str, pickled_key: str, model_case: ModelCase) -> List[str]:
    """loads the validation genes from a pickled file"""
    pickled_genes_path = pickled_genes_path if os.path.isfile(pickled_genes_path) else make_absolute_path("../../vitis_data", pickled_genes_path, start_file=__file__)
    if model_case in [ModelCase.SSR, ModelCase.SSC]:
        with open(pickled_genes_path, 'rb') as handle:
            validation_genes = pickle.load(handle)
            validation_genes = validation_genes[pickled_key]
    else:
        validation_genes = []
    return validation_genes

def load_input_files_training(genome_file_name: str, annotation_file_name: str, tpm_file_name: str, model_case: ModelCase) -> Tuple[Fasta, pd.DataFrame, pd.DataFrame]:
    """loads the input files for training"""
    loaded = load_input_files(genome_file_name=genome_file_name, annotation_file_name=annotation_file_name, tpm_counts_file_name=tpm_file_name, model_case=str(model_case))
    return loaded["genome"], loaded["tpms"], loaded["annotation"]

def set_up_train_val_split_variables(annotation: pd.DataFrame, validation_genes_count: int, validation_fraction: float = 0.2) -> Tuple[int, int]:
    """sets up the variables for the train val split - only used when train_val_split=True"""
    total_sequences = len(annotation)
    target_val_size = int(total_sequences * validation_fraction)  # 20% for validation
    target_train_size = total_sequences - target_val_size  # 80% for training
    
    print(f"=== TRAIN/VAL SPLIT CONFIGURATION ===")
    print(f"Total genes: {total_sequences}")
    print(f"Target validation size: {target_val_size} ({validation_fraction*100:.1f}%)")
    print(f"Target training size: {target_train_size} ({(1-validation_fraction)*100:.1f}%)")
    print("======================================")
    
    return target_val_size, target_train_size

def track_excluded_genes(excluded_genes: List[str], val_chromosome: str) -> None:
    """Track genes consistently excluded across all models"""
    global GLOBAL_EXCLUDED_GENES
    
    # Filter out expected validation chromosome genes and missing TPM genes (these are normal exclusions)
    truly_excluded = [gene for gene in excluded_genes 
                     if not (gene.startswith('chr') or 'missing_tpm' in gene.lower())]
    
    if GLOBAL_EXCLUDED_GENES['model_count'] == 0:
        # First model - initialize with all excluded genes
        GLOBAL_EXCLUDED_GENES['consistently_excluded'] = set(truly_excluded)
    else:
        # Subsequent models - keep only genes excluded by ALL models
        GLOBAL_EXCLUDED_GENES['consistently_excluded'] &= set(truly_excluded)
    
    GLOBAL_EXCLUDED_GENES['model_count'] += 1

def write_final_exclusion_summary(time_stamp: str, output_name: str) -> None:
    """Write summary of genes consistently excluded across ALL models"""
    global GLOBAL_EXCLUDED_GENES
    
    if GLOBAL_EXCLUDED_GENES['model_count'] == 0:
        return
    
    summary_file = make_absolute_path("../../out", "training", f"consistently_excluded_genes_{output_name}_{time_stamp}.txt", start_file=__file__)
    
    with open(summary_file, 'w') as f:
        f.write("=== GENES CONSISTENTLY EXCLUDED FROM ALL MODELS ===\n")
        f.write(f"Training run: {output_name}\n")
        f.write(f"Timestamp: {time_stamp}\n")
        f.write(f"Total models trained: {GLOBAL_EXCLUDED_GENES['model_count']}\n")
        f.write(f"Genes excluded from ALL models: {len(GLOBAL_EXCLUDED_GENES['consistently_excluded'])}\n\n")
        
        if GLOBAL_EXCLUDED_GENES['consistently_excluded']:
            f.write("Consistently excluded genes:\n")
            for gene in sorted(GLOBAL_EXCLUDED_GENES['consistently_excluded']):
                f.write(f"{gene}\n")
        else:
            f.write("No genes were consistently excluded from all models.\n")
    
    print(f"\n=== FINAL EXCLUSION SUMMARY ===")
    print(f"Models trained: {GLOBAL_EXCLUDED_GENES['model_count']}")
    print(f"Consistently excluded genes: {len(GLOBAL_EXCLUDED_GENES['consistently_excluded'])}")
    if len(GLOBAL_EXCLUDED_GENES['consistently_excluded']) > 0:
        print(f"Summary written to: {summary_file}")
    print("=" * 40)

# Keep the chromosome-wise sampling function for potential future use with random splits
def apply_chromosome_wise_sampling(annotation: pd.DataFrame, tpms: pd.DataFrame, validation_genes: List[str], 
                                   val_chromosome: str, validation_fraction: float = 0.2, time_stamp: str = "") -> pd.DataFrame:
    """Apply chromosome-wise sampling - ONLY used when train_val_split=True"""
    print("NOTE: Chromosome-wise sampling enabled for random split mode")
    
    # Store original annotation for comparison
    original_annotation = annotation.copy()
    
    # Pre-filter chromosomes: only keep chromosomes with genes that have TPM data
    valid_chromosomes = []
    invalid_chromosome_stats = {}
    
    for chrom in annotation['Chromosome'].unique():
        chr_genes = annotation[annotation['Chromosome'] == chrom]
        chr_with_tpm = chr_genes.merge(tpms[['target']], left_on='gene_id', right_index=True, how='inner')
        
        if len(chr_with_tpm) > 0:  # Keep chromosomes with any TPM data
            valid_chromosomes.append(chrom)
        else:
            invalid_chromosome_stats[chrom] = len(chr_genes)
            print(f"Excluding chromosome {chrom}: {len(chr_genes)} genes with no TPM data")
    
    # Filter annotation to only valid chromosomes before calculating targets
    annotation = annotation[annotation['Chromosome'].isin(valid_chromosomes)]
    
    # Calculate global targets for 20/80 split (now on filtered data)
    total_genes = len(annotation)
    global_validation_target = len(validation_genes)  # Use actual validation gene count
    global_training_target = int(total_genes * (1 - validation_fraction))  # 80% for training
    
    print(f"=== CHROMOSOME-WISE SAMPLING (RANDOM SPLIT MODE) ===")
    print(f"Valid chromosomes: {len(valid_chromosomes)}")
    print(f"Total genes for sampling: {total_genes}")
    print(f"Global validation target: {global_validation_target}")
    print(f"Global training target: {global_training_target}")
    print("=" * 55)
    
    sampled_genes = []
    
    chromosomes = annotation['Chromosome'].unique()
    for chrom in chromosomes:
        # Get genes for this chromosome
        chr_genes = annotation[annotation['Chromosome'] == chrom]
        total_chr_genes = len(chr_genes)
        
        # Calculate this chromosome's proportional allocation of the global training target
        chr_proportion = total_chr_genes / total_genes
        chr_training_allocation = int(global_training_target * chr_proportion)
        
        # Split into validation and non-validation genes
        chr_validation_genes = chr_genes[chr_genes['gene_id'].isin(validation_genes)]
        chr_non_validation_genes = chr_genes[~chr_genes['gene_id'].isin(validation_genes)]
        
        # Always keep validation genes
        chr_sampled = [chr_validation_genes]
        
        if len(chr_non_validation_genes) > 0:
            # Get TPM targets for non-validation genes
            chr_non_val_with_targets = chr_non_validation_genes.merge(
                tpms[['target']], left_on='gene_id', right_index=True, how='inner'
            )
            
            # Separate by expression level
            high_expr = chr_non_val_with_targets[chr_non_val_with_targets['target'] == 1]  # High expression
            med_expr = chr_non_val_with_targets[chr_non_val_with_targets['target'] == 2]   # Medium expression  
            low_expr = chr_non_val_with_targets[chr_non_val_with_targets['target'] == 0]   # Low expression
            
            # Always keep high and low expression genes
            chr_sampled.extend([high_expr, low_expr])
            
            # Calculate remaining allocation after keeping validation, high, and low expression genes
            current_kept = len(chr_validation_genes) + len(high_expr) + len(low_expr)
            remaining_allocation = chr_training_allocation + len(chr_validation_genes) - current_kept
            
            # Sample from medium expression genes to fill remaining allocation
            if len(med_expr) > 0 and remaining_allocation > 0:
                med_sample_size = min(remaining_allocation, len(med_expr))
                if med_sample_size > 0:
                    sampled_med = med_expr.sample(n=med_sample_size, random_state=42)
                    # Keep only columns that exist in both dataframes
                    common_cols = [col for col in chr_genes.columns if col in sampled_med.columns]
                    chr_sampled.append(sampled_med[common_cols])
        
        # Combine all sampled genes for this chromosome
        if chr_sampled:
            # Only keep original annotation columns for all dataframes
            dfs_to_concat = []
            for df in chr_sampled:
                if not df.empty:
                    # Keep only columns that exist in both df and chr_genes
                    common_cols = [col for col in chr_genes.columns if col in df.columns]
                    dfs_to_concat.append(df[common_cols])
            
            if dfs_to_concat:
                final_chr_sample = pd.concat(dfs_to_concat, ignore_index=True)
            else:
                final_chr_sample = pd.DataFrame(columns=chr_genes.columns)
        else:
            final_chr_sample = pd.DataFrame(columns=chr_genes.columns)
            
        sampled_genes.append(final_chr_sample)
    
    # Check if we have any genes to concatenate
    if not sampled_genes:
        print("ERROR: No genes found after sampling!")
        return pd.DataFrame()
    
    non_empty_genes = [df for df in sampled_genes if not df.empty]
    if not non_empty_genes:
        print("ERROR: All sampled gene dataframes are empty!")
        return pd.DataFrame()
    
    sampled_annotation = pd.concat(non_empty_genes, ignore_index=True)
    
    print(f"Sampling result: {len(original_annotation)} → {len(sampled_annotation)} genes")
    return sampled_annotation

def extract_genes_training(genome_path: str, annotation_path: str, tpm_path: str, extragenic: int, intragenic: int, genes_pickled: str,
                           pickled_key: str, val_chromosome: str, model_case: ModelCase, ignore_small_genes: bool, train_val_split: bool, time_stamp: str,
                           validation_fraction: float, test_specie: str = "") -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """extracts the gene flanking regions for training - simplified for chromosome-based validation"""
    
    if model_case == ModelCase.MSR:
        if test_specie is None:
            raise ValueError("test specie parameter necessary for msr training!")
    
    genome, tpms, annotation = load_input_files_training(genome_path, annotation_path, tpm_path, model_case)
    ssc_training = model_case == ModelCase.SSC        
    validation_genes = set_up_validation_genes(genes_pickled, pickled_key, model_case)
    expected_final_size = 2*(extragenic + intragenic) + 20

    print("=== PRE-FILTERING CHR00 ===")
    original_gene_count = len(annotation)
    valid_chromosomes = []
    
    for chrom in annotation['Chromosome'].unique():
        chr_genes = annotation[annotation['Chromosome'] == chrom]
        chr_with_tpm = chr_genes.merge(tpms[['target']], left_on='gene_id', right_index=True, how='inner')
        
        if len(chr_with_tpm) > 0:  # Keep chromosomes with any TPM data
            valid_chromosomes.append(chrom)
        else:
            print(f"Excluding chromosome {chrom}: {len(chr_genes)} genes with no TPM data")
    
    # Filter annotation to only valid chromosomes
    annotation = annotation[annotation['Chromosome'].isin(valid_chromosomes)]
    filtered_gene_count = len(annotation)

    # For chromosome-based training, NO SAMPLING - use all genes as per original design
    if not train_val_split:
        print(f"=== CHROMOSOME-BASED VALIDATION: chr{val_chromosome} ===")
        print(f"Training genes: ALL genes NOT on chr{val_chromosome}")
        print(f"Validation genes: validation genes ON chr{val_chromosome}")
        print("No gene sampling applied - using original chromosome holdout design")
        print("=" * 60)
    else:
        # Apply sampling only for random splits
        annotation = apply_chromosome_wise_sampling(
            annotation=annotation, 
            tpms=tpms, 
            validation_genes=validation_genes, 
            val_chromosome=val_chromosome, 
            validation_fraction=validation_fraction,
            time_stamp=time_stamp
        )

    # Set up split variables (only affects random splits due to calculate_conditions logic)
    current_val_size, current_train_size = 0, 0
    target_val_size, target_train_size = set_up_train_val_split_variables(
        annotation=annotation, 
        validation_genes_count=len(validation_genes), 
        validation_fraction=validation_fraction
    )
        
    train_seqs, val_seqs, train_targets, val_targets = [], [], [], []
    skipped_genes = []
    
    if len(annotation.columns) == 5:
        annotation["species"] = "specie unknown"
        annotation = annotation[['species', 'Chromosome', 'Start', 'End', 'Strand', 'gene_id']]
        
    for specie, chrom, start, end, strand, gene_id in annotation.values:
        if gene_id not in tpms.index:
            skipped_genes.append(gene_id)
            continue
            
        include_in_validation_set, include_in_training_set = calculate_conditions(
            val_chromosome, model_case, train_val_split, test_specie,
            validation_genes, current_val_size, current_train_size,
            target_val_size, target_train_size, specie, chrom, gene_id
        )

        seq = extract_gene(
            genome=genome, extragenic=extragenic, intragenic=intragenic, ignore_small_genes=ignore_small_genes,
            expected_final_size=expected_final_size, chrom=chrom, start=start, end=end, strand=strand, 
            ssc_training=ssc_training, val_chromosome=val_chromosome
        )
            
        added_val, added_train = append_sequence_training(
            include_as_validation_gene=include_in_validation_set, include_as_training_gene=include_in_training_set, 
            train_targets=train_targets, expected_final_size=expected_final_size, train_seqs=train_seqs, 
            val_seqs=val_seqs, val_targets=val_targets, tpms=tpms, gene_id=gene_id, seq=seq
        )
        
        if added_train == 0 and added_val == 0:
            skipped_genes.append(gene_id)
                
        current_val_size += added_val
        current_train_size += added_train

    # Only check ratio for random splits
    if train_val_split and current_val_size < target_val_size:
        raise ValueError(f"Validation set is not {validation_fraction*100}%. Current size: {current_val_size} genes, Target size: {target_val_size} genes.")

    # Track excluded genes for global analysis
    track_excluded_genes(skipped_genes, val_chromosome)
    
    train_seqs, val_seqs, train_targets, val_targets = np.array(train_seqs), np.array(val_seqs), np.array(train_targets), np.array(val_targets)
    
    # Report final statistics
    actual_train_count = len(train_seqs)
    actual_val_count = len(val_seqs)
    total_used = actual_train_count + actual_val_count
    
    print(f"\n=== FINAL MODEL STATISTICS ===")
    print(f"Training genes: {actual_train_count}")
    print(f"Validation genes: {actual_val_count}")
    print(f"Total genes: {total_used}")
    print(f"Skipped genes: {len(skipped_genes)}")
    print("=" * 35)
    
    if train_seqs.size == 0 or val_seqs.size == 0:
        raise ValueError("Validation sequences or training sequences are empty.")

    mask_sequences(train_seqs=train_seqs, val_seqs=val_seqs, extragenic=extragenic, intragenic=intragenic)
    return train_seqs, train_targets, val_seqs, val_targets

def mask_sequences(train_seqs: np.ndarray, val_seqs: np.ndarray, extragenic: int, intragenic: int):
    """masking the start and end codon of the sequence as zeros, according to original paper"""
    # Masking
    train_seqs[:, extragenic:extragenic + 3, :] = 0
    train_seqs[:, extragenic + (intragenic * 2) + 17:extragenic + (intragenic * 2) + 20, :] = 0
    val_seqs[:, extragenic:extragenic + 3, :] = 0
    val_seqs[:, extragenic + (intragenic * 2) + 17:extragenic + (intragenic * 2) + 20, :] = 0

def deep_cre(x_train: np.ndarray, y_train: np.ndarray, x_val: np.ndarray, y_val: np.ndarray, output_name: str,
             model_case: ModelCase, chrom: str, time_stamp: str, test: bool = False) -> Tuple[float, float, float, float]:
    """trains a deep learning model for CRE prediction"""
    input_seq = Input(shape=(x_train.shape[1], x_train.shape[2]))

    # Conv block 1
    conv = Conv1D(filters=64, kernel_size=8, padding='same')(input_seq)
    conv = Activation('relu')(conv)
    conv = Conv1D(filters=64, kernel_size=8, padding='same')(conv)
    conv = Activation('relu')(conv)
    conv = MaxPool1D(pool_size=8, padding='same')(conv)
    conv = Dropout(0.25)(conv)

    # Conv block 2 and 3
    for n_filters in [128, 64]:
        conv = Conv1D(filters=n_filters, kernel_size=8, padding='same')(conv)
        conv = Activation('relu')(conv)
        conv = Conv1D(filters=n_filters, kernel_size=8, padding='same')(conv)
        conv = Activation('relu')(conv)
        conv = MaxPool1D(pool_size=8, padding='same')(conv)
        conv = Dropout(0.25)(conv)

    # Fully connected block
    output = Flatten()(conv)
    output = Dense(128)(output)
    output = Activation('relu')(output)
    output = Dropout(0.25)(output)
    output = Dense(64)(output)
    output = Activation('relu')(output)
    output = Dense(1)(output)
    output = Activation('sigmoid')(output)

    model = Model(inputs=input_seq, outputs=output)
    model.summary()

    file_name = get_filename_from_path(__file__)
    potential_folder = os.path.dirname(output_name)
    if os.path.exists(potential_folder) and not os.path.isfile(potential_folder):
        checkpoint_path = output_name
    else:
        checkpoint_path = make_absolute_path("../../out/saved_models", f"{output_name}_{chrom}_{model_case}_{file_name}_{time_stamp}.h5", start_file=__file__)
    
    model_chkpt = ModelCheckpoint(filepath=checkpoint_path, save_best_only=True, verbose=1)
    early_stop = EarlyStopping(patience=3) if test else EarlyStopping(patience=10)
    reduce_lr = ReduceLROnPlateau(patience=3, factor=0.1) if test else ReduceLROnPlateau(patience=5, factor=0.1)
    
    model.compile(loss='binary_crossentropy', optimizer=Adam(0.0001),
                  metrics=['accuracy', AUC(curve="ROC", name='auROC'), AUC(curve="PR", name='auPR')])
    model.fit(x_train, y_train, batch_size=64, epochs=100, validation_data=(x_val, y_val),
              callbacks=[early_stop, model_chkpt, reduce_lr])

    loaded_model = load_model(checkpoint_path)
    output = loaded_model.evaluate(x_val, y_val)
    return output

def balance_dataset(x: np.ndarray, y: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """randomly down samples the majority class to balance the dataset"""
    # Random down sampling to balance data
    low_train, high_train = np.where(y == 0)[0], np.where(y == 1)[0]
    min_class = min([len(low_train), len(high_train)])
    selected_low_train = np.random.choice(low_train, min_class, replace=False)
    selected_high_train = np.random.choice(high_train, min_class, replace=False)
    x_train = np.concatenate([
        np.take(x, selected_low_train, axis=0),
        np.take(x, selected_high_train, axis=0)
    ], axis=0)
    y_train = np.concatenate([
        np.take(y, selected_low_train, axis=0),
        np.take(y, selected_high_train, axis=0)
    ], axis=0)
    x_train, y_train = shuffle(x_train, y_train, random_state=42)
    return x_train, y_train

def train_deep_cre(genome_path: str, annotation_path: str, tpm_path: str, extragenic: int, intragenic: int, genes_picked: str,
                   val_chromosome: str, output_name: str, model_case: ModelCase, ignore_small_genes: bool, train_val_split: bool, time_stamp: str,
                   validation_fraction: float,  test_specie: Optional[str] = None, pickled_key: Optional[str] = None, test: bool = False) -> Tuple[float, float, float, float]:
    """trains the deepCRE model"""
    train_seqs, train_targets, val_seqs, val_targets = extract_genes_training(
        genome_path, annotation_path, tpm_path, extragenic, intragenic,
        genes_picked, pickled_key, val_chromosome, model_case, ignore_small_genes,
        train_val_split=train_val_split, test_specie=test_specie, time_stamp=time_stamp, validation_fraction=validation_fraction
    )
    
    x_train, y_train = balance_dataset(train_seqs, train_targets)
    x_val, y_val = balance_dataset(val_seqs, val_targets)
    
    output = deep_cre(x_train=x_train, y_train=y_train, x_val=x_val, y_val=y_val, output_name=output_name,
                      model_case=model_case, chrom=val_chromosome, time_stamp=time_stamp, test=test)
    return output

def get_chromosomes(species_info: List[Dict[str, Any]]) -> List[str]:
    """loads the chromosomes from the species info"""
    chromosomes = species_info[0]["chromosomes"]
    if not chromosomes:
        raise ValueError("For SSR/SSC training, a list of chromosomes must be provided!")
    return chromosomes

def run_ssr(species_info: List[Dict[str, Any]], general_info: Dict[str, Any], time_stamp: str, test: bool = False) -> List[Dict[str, Any]]:
    """Runs model training for a single species training run with automatic resume"""
    specie_info = species_info[0]
    print(f"Single species Training on genome:\n------------------------------")
    print(specie_info["genome"])
    print("------------------------------\n")

    all_chromosomes = get_chromosomes(species_info)

    # Automatically filter out already completed chromosomes
    chromosomes = filter_remaining_chromosomes(
        all_chromosomes, general_info["output_name"], general_info["model_case"]
    )
    
    if not chromosomes:
        print("All chromosomes already trained! No work to do.")
        return []
    
    print(f"Chromosomes to train: {chromosomes}")
    combined_results = []

    for val_chrom in chromosomes:
        print(f"Using chromosome {val_chrom} as validation chromosome")
        
        try:
            results = train_deep_cre(
                genome_path=specie_info["genome"], annotation_path=specie_info["annotation"], tpm_path=specie_info["targets"],
                extragenic=general_info["extragenic"], intragenic=general_info["intragenic"], genes_picked=specie_info["pickle_file"],
                val_chromosome=val_chrom, output_name=general_info["output_name"], model_case=general_info["model_case"],
                pickled_key=specie_info["pickle_key"], ignore_small_genes=general_info["ignore_small_genes"],
                train_val_split=general_info["train_val_split"], time_stamp=time_stamp, test=test, validation_fraction=general_info["validation_fraction"]
            )

            results_with_info = {
                "loss": results[0], "accuracy": results[1], "auROC": results[2], "auPR": results[3], "test": val_chrom,
            }
            combined_results.append(results_with_info)
            print(f"Results for chromosome {val_chrom}: {results}")
            
        except Exception as e:
            print(f"Error training chromosome {val_chrom}: {e}")
            continue
        finally:
            cleanup_memory()
            
    return combined_results

def run_msr(species_info: List[Dict[str, Any]], general_info: Dict[str, Any], time_stamp: str, test: bool = False) -> List[Dict[str, Any]]:
    """runs model training for a multi-species training run"""
    species: List[str] = [specie["species_name"] for specie in species_info]
    naming = "_".join([specie.replace(" ", "").replace(".", "") for specie in species])

    # generate concat files
    tpm_path = combine_tpms(species_data=species_info, naming=naming)
    genome_path = combine_fasta(species_data=species_info, naming=naming)
    annotation_path = combine_annotations(species_data=species_info, naming=naming)
    combined_results = []

    for test_specie_info in species_info:
        train_species = [specie for specie in species if specie != test_specie_info["species_name"]]
        print(f'Training on species: {train_species}')
        print(f'Testing on specie: {test_specie_info["species_name"]}')

        results = train_deep_cre(
            genome_path=genome_path, annotation_path=annotation_path, tpm_path=tmp_path, extragenic=general_info["extragenic"],
            intragenic=general_info["intragenic"], genes_picked=test_specie_info["pickle_file"], val_chromosome=test_specie_info["species_name"],
            output_name=general_info["output_name"], model_case=general_info["model_case"], ignore_small_genes=general_info["ignore_small_genes"],
            train_val_split=general_info["train_val_split"], test_specie=test_specie_info["species_name"], time_stamp=time_stamp, test=test, validation_fraction=0.2
        ) 
        results_with_info = {
            'loss': results[0], 'accuracy': results[1], 'auROC': results[2], 'auPR': results[3], 'test': test_specie_info['species_name'],
        }
        combined_results.append(results_with_info)
        print(f"Results for validation species {test_specie_info['species_name']}: {results}")
        cleanup_memory()
    return combined_results

def parse_args() -> argparse.Namespace:
    """parses the arguments provided by the user"""
    parser = argparse.ArgumentParser(
        prog='deepCRE',
        description="This script performs deepCRE training."
    )
    parser.add_argument('--input', "-i", required=True,
                        help="json file containing the required input parameters")
    return parser.parse_args()

def parse_input_file(input_file: str) -> Tuple[ParsedInputs, List[Tuple[str, int, Exception]], int]:
    """parses the input file"""
    possible_general_parameters = {
        "model_case": None, "genome": None, "annotation": None, "targets": None, "output_name": None,
        "chromosomes": "", "pickle_key": None, "pickle_file": "validation_genes.pickle",
        "ignore_small_genes": True, "train_val_split": False, "extragenic": 1000, "intragenic": 500,
        "validation_fraction": 0.2,
    }

    possible_species_parameters = {
        "genome": None, "annotation": None, "targets": None, "chromosomes": "",
        "pickle_key": None, "pickle_file": "validation_genes.pickle", "species_name": None,
    }
    
    inputs, failed_trainings, input_length = ParsedInputs.parse(
        input_file, possible_general_parameters=possible_general_parameters, 
        possible_species_parameters=possible_species_parameters, multiple_species_required_msr=True
    )
    inputs = inputs.replace_both()
    return inputs, failed_trainings, input_length

def train_models(inputs: ParsedInputs, failed_trainings: List[Tuple[str, int, Exception]], input_length: int, test: bool = False) -> List[Tuple[str, int, Exception]]:
    """trains the models for the given inputs"""
    os.environ['TF_FORCE_GPU_ALLOW_GROWTH'] = 'true'
    training_results_path = make_absolute_path('../../out', "training", start_file=__file__)
    models_path = make_absolute_path("../../out/saved_models", start_file=__file__)
    os.makedirs(training_results_path, exist_ok=True)
    os.makedirs(models_path, exist_ok=True)

    file_name = get_filename_from_path(__file__)
    time_stamp = get_time_stamp()
    
    # Initialize global tracking
    global GLOBAL_EXCLUDED_GENES
    GLOBAL_EXCLUDED_GENES['timestamp'] = time_stamp
    
    for i, run_info in enumerate(inputs):
        try:
            gen_info = run_info.general_info
            spec_info = run_info.species_info
            
            if run_info.is_msr():
                results = run_msr(species_info=spec_info, general_info=gen_info, time_stamp=time_stamp, test=test)
            else:
                results = run_ssr(species_info=spec_info, general_info=gen_info, time_stamp=time_stamp, test=test)
            
            if results:  # Only save if we have results
                results_df = pd.DataFrame(results, columns=['test', 'loss', 'accuracy', 'auROC', 'auPR'])
                save_file = os.path.join(training_results_path, f"{gen_info['output_name']}_{file_name}_{gen_info['model_case']}_{time_stamp}.csv")
                results_df.to_csv(path_or_buf=save_file, index=False)
                print(results_df.head())
            
            # Write final exclusion summary
            write_final_exclusion_summary(time_stamp, gen_info['output_name'])
            
        except TerminationError as e:
            print(run_info)
            raise e
        except Exception as e:
            print(f"Error in training: {e}")
            print(run_info)
            failed_trainings.append((run_info.general_info["output_name"], i, e))
    
    result_summary(failed_runs=failed_trainings, input_length=input_length, script=get_filename_from_path(__file__))
    return failed_trainings

def main() -> None:
    """main function for training the models"""
    configure_gpu_memory()
    args = parse_args()
    inputs, failed_trainings, input_length = parse_input_file(args.input)
    train_models(inputs, failed_trainings, input_length)

if __name__ == "__main__":
    main()