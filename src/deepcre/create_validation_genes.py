import os
import pickle
from typing import Tuple
import pandas as pd
from Bio import SeqIO
import re
import argparse
from vitis_cre.src.deepcre.utils import make_absolute_path


def parse_args() -> argparse.Namespace:
    """Parse command line arguments.

    Returns:
        argparse.Namespace: parsed arguments
    """
    parser = argparse.ArgumentParser(
                        prog='create validation genes',
                        description="This script can be used to create a list of validation genes for a given proteome.")
    parser.add_argument('--proteins', "-p", help="Path to Proteome file", required=True)
    parser.add_argument('--blast_outputs', "-b", help="Path to BLAST output file. Needs to be in outfmt 6", required=True)
    parser.add_argument('--pickle_key_name', "-k", help="Key for the pickle dictionary", required=True)
    parser.add_argument('--output', "-o", help="Path to an output file. By default will create a file in the src/deepCRE folder.", required=False, default="")

    args = parser.parse_args()
    return args


def get_save_path(output_path: str, pickle_key: str) -> str:
    """Get the path to save the pickle file.

    creates path to generated name based on the pickle_key. If output_path is provided, it will be used as the output path.

    Args:
        output_path (str): optional user provided output path
        pickle_key (str): key for the pickle dictionary

    Raises:
        FileExistsError: If the file already exists at the output path

    Returns:
        str: path to save the pickle file
    """
    if not output_path:
        output_path = make_absolute_path("../../vitis_data", f"validation_genes_{pickle_key}.pickle", start_file=__file__)

    # Remove any existing pickle file
    if os.path.exists(output_path):
        raise FileExistsError(f"File already exists at {output_path}. Please remove it or provide a different output path.")
    return output_path


def get_queries(description: str) -> Tuple[str, str]:
    """Get the queries to extract gene and chromosome information from the protein description.

    Args:
        description (str): Description of the protein

    Returns:
        Tuple[str, str]: name and separator for gene and chromosome information
    """
    if "gene=" in description:
        gene_query = "gene="
    elif "gene:" in description:
        gene_query = "gene:"
    else:
        gene_query = ""
    
    if "seq_id=" in description:
        chrom_query = "seq_id="
    elif "chromosome:" in description:
        chrom_query = "chromosome:"
    else:
        chrom_query = ""
    return gene_query, chrom_query


def main() -> None:
    """Main function for creating validation genes.

    This function reads the proteome file and BLAST output file to create a list of genes that have no
    homologs in other chromosomes. The list is saved as an entry in a dictionary under the provided key.
    """
    args = parse_args()
    proteome_path = args.proteins
    blast_output_path = args.blast_outputs
    pickle_key_name = args.pickle_key_name # Key for the pickle dictionary
    output_path = args.output
    print(proteome_path, blast_output_path, pickle_key_name, output_path)

    pickle_file_path = get_save_path(output_path=output_path, pickle_key=pickle_key_name)

    # Dictionary to store validation genes
    validation_genes = {}

    if os.path.exists(blast_output_path):
        # Parse protein sequences and build a DataFrame
        info = []
        for rec in SeqIO.parse(proteome_path, 'fasta'):
            description = rec.description
            gene_query, chrom_query = get_queries(description)
            if not (gene_query and chrom_query):
                continue
            protein_id = rec.id.split('|')[0]  # Extract protein_id
            gene_match = re.search(f'{gene_query}([^| ]+)', description)
            chrom_match = re.search(f'{chrom_query}([^| ]+)', description)
            gene_id = gene_match.group(1) if gene_match else None
            chrom = chrom_match.group(1) if chrom_match else None

            if chrom is not None and chrom.count(":") == 4:
                chrom = chrom.split(":")[1]
            info.append([gene_id, protein_id, chrom])

    
        info = pd.DataFrame(info, columns=['gene_id', 'protein_id', 'chrom'])
        info.index = info.protein_id.tolist()  #type:ignore

        ### DEBUG CHECKPOINT: Show how many proteins parsed
        print(f"Parsed {len(info)} protein records from the proteome.")
        print("Example protein IDs in index:", info.index[:10].tolist()) #type:ignore


        # Read BLAST output
        blast_out = pd.read_csv(
            blast_output_path,
            sep='\t',
            names=[
                'qseqid', 'sseqid', 'pident', 'length', 'mismatch',
                'gap_open', 'qstart', 'qend', 'sstart', 'send', 'evalue',
                'bitscore'
            ]
        )

        ### DEBUG CHECKPOINT: Show raw BLAST rows
        print(f"Read {len(blast_out)} rows from BLAST output.")
        print("Example qseqids in BLAST:", blast_out['qseqid'].unique()[:10].tolist())

        # Filter BLAST output for valid protein IDs
        valid_ids = [seq_id for seq_id in blast_out.qseqid.tolist() if seq_id in info.index]
        blast_out = blast_out[blast_out['qseqid'].isin(valid_ids)]

        if not blast_out.empty:
            # Map gene and chromosome information
            blast_out['qgene'] = info.loc[blast_out.qseqid.tolist(), 'gene_id'].values
            blast_out['sgene'] = info.loc[blast_out.sseqid.tolist(), 'gene_id'].values
            blast_out['qchrom'] = info.loc[blast_out.qseqid.tolist(), 'chrom'].values
            blast_out['schrom'] = info.loc[blast_out.sseqid.tolist(), 'chrom'].values


            ### DEBUG CHECKPOINT: Check mapping success
            print("Sample mappings:")
            print(blast_out[['qseqid', 'sseqid', 'qgene', 'sgene', 'qchrom', 'schrom']].head())

            # Apply filtering criteria
            pre_filter_len = len(blast_out)
            blast_out = blast_out[blast_out['evalue'] < 0.001]
            blast_out = blast_out[blast_out['bitscore'] >= 50]
            blast_out = blast_out[~blast_out['schrom'].isin(['Pt', 'Mt'])]
            blast_out = blast_out[~blast_out['qchrom'].isin(['Pt', 'Mt'])]

            ### DEBUG CHECKPOINT: After applying filters
            print(f"BLAST filtered from {pre_filter_len} to {len(blast_out)} rows.")

            # Identify validation genes
            val_set = []
            for gene_name, gene_grp in blast_out.groupby('qgene'):
                query_chrom = gene_grp['qchrom'].iloc[0] 
                if (gene_grp['schrom'] == query_chrom).all():
                    # If all subject chrom match the query chrom, add to validation set
                    # if at least one homolog is on a different chromosome, skip this gene
                    val_set.append(gene_name)

            # Save flat list of validation genes for this output
            validation_genes[pickle_key_name] = val_set
            
            ### DEBUG CHECKPOINT: Show how many validation genes found
            print(f"Identified {len(val_set)} validation genes.")
        else:
            print("No BLAST matches remained after filtering based on valid protein IDs")

        # Count how many validation genes per chromosome
        val_chr_counts = info[info['gene_id'].isin(val_set)]['chrom'].value_counts().sort_index()

        print("\nValidation gene counts per chromosome:")
        for chrom, count in val_chr_counts.items():
            print(f"{chrom}: {count}")



    else:
        print(f"BLAST output file not found at {blast_output_path}")

    # Save the dictionary with the flattened list of gene IDs
    with open(pickle_file_path, 'wb') as pickle_file:
        pickle.dump(validation_genes, pickle_file, protocol=pickle.HIGHEST_PROTOCOL)

    print(f"Unique chromosomes in val set: {info[info['gene_id'].isin(val_set)]['chrom'].unique()}")
    print(f"Validation genes saved to {pickle_file_path}")

if __name__ == "__main__":
    main()