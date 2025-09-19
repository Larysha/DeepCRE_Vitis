import argparse
import os
from typing import List, Tuple, Dict, Any
import pandas as pd
import numpy as np
from pyfaidx import Fasta

from vitis_cre.src.deepcre.utils import make_absolute_path, load_input_files, get_filename_from_path, get_time_stamp
from vitis_cre.src.deepcre.parsing import ParsedInputs, RunInfo

"""
  Usage Examples:

  Single analysis:
  python create_generic_features.py --genome vitis_vinifera_PN40024_dna.fa \
    --annotation vitis_vinifera_PN40024.gff3 --output vitis_analysis

  Batch processing:
  python -m vitis_cre.src.create_generic_features \
    --input vitis_cre/src/inputs/vitis_training.json
"""


def gc_content(sequence: str) -> float:
    """Calculate GC content of a DNA sequence.

    Args:
        sequence (str): DNA sequence

    Returns:
        float: GC content as a proportion (0-1)
    """
    sequence = str(sequence).upper()
    gc_count = sequence.count('G') + sequence.count('C')
    return gc_count / len(sequence) if len(sequence) > 0 else 0.0


def cpg_percentage(sequence: str) -> float:
    """Calculate CpG dinucleotide percentage in a DNA sequence.

    Args:
        sequence (str): DNA sequence

    Returns:
        float: CpG percentage (0-100)
    """
    sequence = str(sequence).upper()
    cpg_count = sequence.count('CG')
    return (cpg_count / len(sequence)) * 100 if len(sequence) > 0 else 0.0


def get_proximal_promoter_terminator_features(genome: Fasta, annotation: pd.DataFrame, 
                                              target_chromosomes: Tuple[str, ...] = (), 
                                              flank_length: int = 1000) -> pd.DataFrame:
    """Extract GC content and CpG percentage from promoter and terminator regions.

    Args:
        genome (Fasta): Genome sequence data
        annotation (pd.DataFrame): Gene annotation data
        target_chromosomes (Tuple[str, ...]): Chromosomes to analyze (empty = all)
        flank_length (int): Length of flanking regions to extract

    Returns:
        pd.DataFrame: DataFrame with promoter and terminator features
    """
    features = []
    
    for _, row in annotation.iterrows():
        chrom = row['Chromosome']
        start = row['Start']
        end = row['End']
        strand = row['Strand']
        gene_id = row['gene_id']
        
        # Skip if not in target chromosomes (empty tuple means all chromosomes)
        if target_chromosomes and chrom not in target_chromosomes:
            continue
            
        # Skip if chromosome not in genome
        if chrom not in genome.keys():
            continue
        
        try:
            if strand == '-':
                # For negative strand: promoter is downstream of gene end
                prom_start = max(0, end)
                prom_end = min(len(genome[chrom]), end + flank_length)
                # Terminator is upstream of gene start
                term_start = max(0, start - flank_length)
                term_end = min(len(genome[chrom]), start)
            else:
                # For positive strand: promoter is upstream of gene start
                prom_start = max(0, start - flank_length)
                prom_end = min(len(genome[chrom]), start)
                # Terminator is downstream of gene end
                term_start = max(0, end)
                term_end = min(len(genome[chrom]), end + flank_length)
            
            # Extract sequences
            promoter_seq = genome[chrom][prom_start:prom_end]
            terminator_seq = genome[chrom][term_start:term_end]
            
            # For negative strand, get reverse complement
            if strand == '-':
                promoter_seq = promoter_seq.reverse.complement
                terminator_seq = terminator_seq.reverse.complement
            
            # Calculate features
            feature_row = {
                'gene_id': gene_id,
                'Chromosome': chrom,
                'GC_promoter': gc_content(str(promoter_seq)),
                'CpG_promoter': cpg_percentage(str(promoter_seq)),
                'GC_terminator': gc_content(str(terminator_seq)),
                'CpG_terminator': cpg_percentage(str(terminator_seq)),
                'promoter_length': len(str(promoter_seq)),
                'terminator_length': len(str(terminator_seq))
            }
            features.append(feature_row)
            
        except Exception as e:
            print(f"Warning: Could not extract features for gene {gene_id} on {chrom}: {e}")
            continue
    
    return pd.DataFrame(features)


def get_utr_features(genome: Fasta, annotation_path: str, target_chromosomes: Tuple[str, ...] = ()) -> pd.DataFrame:
    """Extract UTR length, GC content, and CpG percentage features.
    
    Note: This function requires UTR annotations in the GFF/GTF file.
    Many genome annotations may not include detailed UTR information.

    Args:
        genome (Fasta): Genome sequence data
        annotation_path (str): Path to annotation file with UTR features
        target_chromosomes (Tuple[str, ...]): Chromosomes to analyze (empty = all)

    Returns:
        pd.DataFrame: DataFrame with UTR features
    """
    try:
        # Load annotation with UTR features
        import pyranges as pr
        if annotation_path.endswith('.gtf'):
            gene_models = pr.read_gtf(annotation_path, as_df=True)
        else:
            gene_models = pr.read_gff3(annotation_path, as_df=True)
        
        # Filter for UTR features
        utr_models = gene_models[gene_models['Feature'].isin(['five_prime_utr', 'three_prime_utr'])]
        
        if utr_models.empty:
            print("Warning: No UTR annotations found in the file. Returning empty DataFrame.")
            return pd.DataFrame()
        
        # Filter by chromosomes if specified
        if target_chromosomes:
            utr_models = utr_models[utr_models['Chromosome'].isin(target_chromosomes)]
        
        # Get unique genes
        all_genes = utr_models['gene_id'].unique() if 'gene_id' in utr_models.columns else []
        
        features = []
        
        for gene_id in all_genes:
            gene_utrs = utr_models[utr_models['gene_id'] == gene_id]
            
            # Initialize feature dictionary
            feature_row = {
                'gene_id': gene_id,
                "5_UTR_length": 0,
                "GC_5_UTR": 0.0,
                "CpG_5_UTR": 0.0,
                "3_UTR_length": 0,
                "GC_3_UTR": 0.0,
                "CpG_3_UTR": 0.0
            }
            
            # Process each UTR type
            for utr_type in ['five_prime_utr', 'three_prime_utr']:
                utr_data = gene_utrs[gene_utrs['Feature'] == utr_type]
                
                if not utr_data.empty:
                    # Take the longest UTR if multiple exist
                    utr_data = utr_data.copy()
                    utr_data['length'] = utr_data['End'] - utr_data['Start']
                    longest_utr = utr_data.loc[utr_data['length'].idxmax()]
                    
                    chrom = longest_utr['Chromosome']
                    start = longest_utr['Start']
                    end = longest_utr['End']
                    strand = longest_utr['Strand']
                    
                    try:
                        if chrom in genome.keys():
                            utr_seq = genome[chrom][start:end]
                            if strand == '-':
                                utr_seq = utr_seq.reverse.complement
                            
                            if utr_type == 'five_prime_utr':
                                feature_row["5_UTR_length"] = end - start
                                feature_row["GC_5_UTR"] = gc_content(str(utr_seq))
                                feature_row["CpG_5_UTR"] = cpg_percentage(str(utr_seq))
                            else:  # three_prime_utr
                                feature_row["3_UTR_length"] = end - start
                                feature_row["GC_3_UTR"] = gc_content(str(utr_seq))
                                feature_row["CpG_3_UTR"] = cpg_percentage(str(utr_seq))
                    except Exception as e:
                        print(f"Warning: Could not extract {utr_type} for gene {gene_id}: {e}")
            
            features.append(feature_row)
        
        return pd.DataFrame(features)
        
    except Exception as e:
        print(f"Error processing UTR features: {e}")
        return pd.DataFrame()


def generate_genomic_features(genome_file_name: str, annotation_file_name: str, 
                             target_chromosomes: Tuple[str, ...] = (), 
                             flank_length: int = 1000, 
                             output_name: str = "") -> pd.DataFrame:
    """Generate comprehensive genomic features for genes.

    Args:
        genome_file_name (str): Name or path of genome FASTA file
        annotation_file_name (str): Name or path of annotation file
        target_chromosomes (Tuple[str, ...]): Chromosomes to analyze (empty = all)
        flank_length (int): Length of flanking regions for promoter/terminator analysis
        output_name (str): Base name for output files

    Returns:
        pd.DataFrame: Combined genomic features
    """
    # Load input files
    loaded_files = load_input_files(
        genome_file_name=genome_file_name,
        annotation_file_name=annotation_file_name
    )
    
    genome = loaded_files["genome"]
    annotation = loaded_files["annotation"]
    
    print(f"Generating genomic features for {len(annotation)} genes...")
    
    # Get promoter and terminator features
    print("Extracting promoter and terminator features...")
    prox_features = get_proximal_promoter_terminator_features(
        genome=genome,
        annotation=annotation,
        target_chromosomes=target_chromosomes,
        flank_length=flank_length
    )
    
    # Get UTR features (if available in annotation)
    print("Extracting UTR features...")
    annotation_path = annotation_file_name if os.path.isfile(annotation_file_name) else make_absolute_path("../../vitis_data/gene_models", annotation_file_name, start_file=__file__)
    utr_features = get_utr_features(
        genome=genome,
        annotation_path=annotation_path,
        target_chromosomes=target_chromosomes
    )
    
    # Combine features
    if not utr_features.empty:
        combined_features = prox_features.merge(utr_features, on='gene_id', how='left')
        print(f"Combined {len(prox_features)} promoter/terminator features with {len(utr_features)} UTR features")
    else:
        combined_features = prox_features
        print("No UTR features found - using only promoter/terminator features")
    
    # Save results
    if output_name:
        output_dir = make_absolute_path("../../out", "genomic_features", start_file=__file__)
        os.makedirs(output_dir, exist_ok=True)
        
        timestamp = get_time_stamp()
        output_file = os.path.join(output_dir, f"{output_name}_genomic_features_{timestamp}.csv")
        combined_features.to_csv(output_file, index=False)
        print(f"Genomic features saved to: {output_file}")
    
    return combined_features


def run_feature_generation(inputs: ParsedInputs, failed_runs: List[Tuple], input_length: int) -> List[Tuple]:
    """Run genomic feature generation for multiple inputs.

    Args:
        inputs (ParsedInputs): Parsed input configurations
        failed_runs (List[Tuple]): List of failed runs
        input_length (int): Total number of input runs

    Returns:
        List[Tuple]: Updated list of failed runs
    """
    for i, run_info in enumerate(inputs):
        try:
            output_name = run_info.general_info.get("output_name", f"run_{i}")
            target_chromosomes = tuple(run_info.general_info.get("chromosomes", []))
            flank_length = run_info.general_info.get("extragenic", 1000)
            
            print(f"\n--- Processing run {i+1}/{input_length}: {output_name} ---")
            
            features = generate_genomic_features(
                genome_file_name=run_info.general_info["genome"],
                annotation_file_name=run_info.general_info["annotation"],
                target_chromosomes=target_chromosomes,
                flank_length=flank_length,
                output_name=output_name
            )
            
            print(f"Generated {len(features)} feature rows for {output_name}")
            
        except Exception as e:
            print(f"Error processing run {i+1}: {e}")
            failed_runs.append((output_name, i, e))
    
    return failed_runs


def parse_input_file(file_path: str) -> Tuple[ParsedInputs, List[Tuple], int]:
    """Parse input file for feature generation.

    Args:
        file_path (str): Path to input JSON file

    Returns:
        Tuple[ParsedInputs, List[Tuple], int]: Parsed inputs, failed runs, input length
    """
    possible_general_parameters = {
        "genome": None,
        "annotation": None,
        "output_name": None,
        "chromosomes": [],
        "extragenic": 1000,
    }
    
    possible_species_parameters = {
        "chromosomes": [],
    }
    
    inputs, failed_runs, input_length = ParsedInputs.parse(
        file_path, 
        possible_general_parameters=possible_general_parameters,
        possible_species_parameters=possible_species_parameters
    )
    inputs = inputs.replace_both()
    return inputs, failed_runs, input_length


def parse_args() -> argparse.Namespace:
    """Parse command line arguments.

    Returns:
        argparse.Namespace: Parsed arguments
    """
    parser = argparse.ArgumentParser(
        prog='create_generic_features',
        description="Generate genomic features (GC content, CpG percentage, UTR lengths) for genes."
    )
    
    parser.add_argument(
        '--input', '-i',
        help="JSON file containing input parameters for feature generation",
        required=True
    )
    
    parser.add_argument(
        '--genome', '-g',
        help="Genome FASTA file (overrides JSON input)",
        required=False
    )
    
    parser.add_argument(
        '--annotation', '-a',
        help="Annotation GFF3/GTF file (overrides JSON input)",
        required=False
    )
    
    parser.add_argument(
        '--output', '-o',
        help="Output name for generated features",
        required=False,
        default=""
    )
    
    parser.add_argument(
        '--flank_length', '-f',
        help="Length of flanking regions for promoter/terminator analysis",
        type=int,
        default=1000
    )
    
    return parser.parse_args()


def main() -> None:
    """Main function for genomic feature generation."""
    args = parse_args()
    
    # If individual files specified, run single analysis
    if args.genome and args.annotation:
        print("Running single genomic feature analysis...")
        features = generate_genomic_features(
            genome_file_name=args.genome,
            annotation_file_name=args.annotation,
            flank_length=args.flank_length,
            output_name=args.output or "single_analysis"
        )
        print(f"Analysis complete. Generated {len(features)} feature rows.")
        print("\nSample features:")
        print(features.head())
        
    else:
        # Parse input file and run multiple analyses
        print("Running batch genomic feature analysis from input file...")
        inputs, failed_runs, input_length = parse_input_file(args.input)
        failed_runs = run_feature_generation(inputs, failed_runs, input_length)
        
        # Summary
        successful_runs = input_length - len(failed_runs)
        print(f"\n--- Summary ---")
        print(f"Successful runs: {successful_runs}/{input_length}")
        
        if failed_runs:
            print("Failed runs:")
            for output_name, run_idx, error in failed_runs:
                print(f"  - {output_name} (run {run_idx + 1}): {error}")


if __name__ == "__main__":
    main()