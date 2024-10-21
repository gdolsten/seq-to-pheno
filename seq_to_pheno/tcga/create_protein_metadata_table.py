import os
import pandas as pd
import argparse
import time
import requests
import numpy as np  # Make sure to import numpy if not already imported


def read_fasta_sequence(fasta_path):
    """
    Reads a FASTA file and returns the sequence as a string.
    Assumes a single sequence per file.
    """
    sequence = ''
    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                continue  # Skip header lines
            sequence += line.strip()
    return sequence

def fetch_gene_names(transcript_ids):
    transcript_ids = list(transcript_ids)
    server = "https://grch37.rest.ensembl.org"
    endpoint = "/lookup/id"
    headers = {"Content-Type": "application/json"}
    gene_name_map = {}

    batch_size = 200  # Reduce batch size if necessary
    for i in range(0, len(transcript_ids), batch_size):
        batch = transcript_ids[i:i+batch_size]
        data = {"ids": batch}
        try:
            response = requests.post(f"{server}{endpoint}", headers=headers, json=data)
            if not response.ok:
                response.raise_for_status()
            decoded = response.json()
            for tid, info in decoded.items():
                if info:
                    gene_name = info.get('display_name')
                    gene_name_map[tid] = gene_name
                else:
                    gene_name_map[tid] = None
        except requests.exceptions.RequestException as e:
            print(f"Error fetching gene names for batch starting at index {i}: {e}")
            # You may choose to retry or skip this batch
            continue
        time.sleep(1)  # To respect rate limits
    return gene_name_map


def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Create protein sequences metadata table.')
    parser.add_argument('--include-sequences', action='store_true',
                        help='Include actual protein sequences in the output instead of file paths.')
    args = parser.parse_args()

    # Paths to directories
    mutated_proteins_dir = 'seq_to_pheno/tcga/data/variants/mutated_proteins'
    wildtype_proteins_dir = 'seq_to_pheno/tcga/data/variants/wildtype_proteins'
    metadata_file = 'seq_to_pheno/tcga/data/new_metadata.csv'
    output_file = 'seq_to_pheno/tcga/data/protein_sequences_metadata.tsv'

    # Read the metadata file into a DataFrame
    metadata_df = pd.read_csv(metadata_file, sep=',', low_memory=False)

    # List all mutated protein FASTA files
    mutated_files = [f for f in os.listdir(mutated_proteins_dir) if f.endswith('.fasta')]

    data_rows = []

    for filename in mutated_files:
        # Extract sample ID and transcript ID from filename
        # Filename format: {sample_id}_{transcript_id}_mutated.fasta
        parts = filename.split('_')
        if len(parts) < 3:
            print(f"Filename {filename} does not match expected format.")
            continue

        sample_id = parts[0]
        transcript_id = parts[1]
        
        # Paths to the mutated and wildtype protein files
        mutated_protein_path = os.path.join(mutated_proteins_dir, filename)
        wildtype_protein_filename = f"{transcript_id}.fasta"
        wildtype_protein_path = os.path.join(wildtype_proteins_dir, wildtype_protein_filename)
        
        # Check if the wildtype protein file exists
        if not os.path.exists(wildtype_protein_path):
            print(f"Wildtype protein file {wildtype_protein_filename} not found.")
            continue

        # Get metadata for the sample
        sample_metadata = metadata_df[metadata_df['wgs_aliquot_id'] == sample_id]

        if sample_metadata.empty:
            print(f"No metadata found for sample {sample_id}.")
            continue

        # Extract required metadata fields
        sample_metadata = sample_metadata.iloc[0]  # Get the first matching row

        # Columns to extract
        columns_to_extract = {
            'aliquot_id': 'aliquot_id',
            'wgs_aliquot_id': 'wgs_aliquot_id',
            'Cancer Type': 'study',
            'Cancer Stage': 'tumour_stage',
            'Donor Survival Time': 'donor_survival_time',
            'Donor Vital Status': 'donor_vital_status',
            'Donor Age at Diagnosis': 'donor_age_at_diagnosis',
            'Tumour Grade': 'tumour_grade',
            'Donor Sex': 'donor_sex',
            'Histology Abbreviation': 'histology_abbreviation',
            # Add more fields as needed
        }

        # Build a dictionary for the row
        row = {
            'aliquot_id': sample_id,
            'transcript_id': transcript_id
        }

        # Include sequences or file paths based on the argument
        if args.include_sequences:
            # Read the sequences from the FASTA files
            mutated_sequence = read_fasta_sequence(mutated_protein_path)
            wildtype_sequence = read_fasta_sequence(wildtype_protein_path)
            row['mutated_protein'] = mutated_sequence
            row['wildtype_protein'] = wildtype_sequence
        else:
            # Include the file paths
            row['mutated_protein'] = mutated_protein_path
            row['wildtype_protein'] = wildtype_protein_path

        for col_name, meta_col in columns_to_extract.items():
            value = sample_metadata.get(meta_col, '')
            row[col_name] = value

        data_rows.append(row)

        # Limit to first 50 entries if including sequences
        # if args.include_sequences and len(data_rows) >= 50:
            # break

    # Create a DataFrame from the collected rows
    result_df = pd.DataFrame(data_rows)
    
    # Ensure transcript IDs are strings
    result_df['transcript_id'] = result_df['transcript_id'].astype(str)
    
    # Clean transcript IDs
    result_df['transcript_id'] = result_df['transcript_id'].str.split('.').str[0].str.strip()

    # Fetch gene names
    transcript_ids = result_df['transcript_id'].unique().tolist()
    gene_name_map = fetch_gene_names(transcript_ids)
    mapping_df = pd.DataFrame(list(gene_name_map.items()), columns=['transcript_id', 'gene_name'])
    metadata_with_gene = result_df.merge(mapping_df, on='transcript_id', how='left')

    # Save the DataFrame to a TSV file
    metadata_with_gene.to_csv(output_file, index=False, sep='\t')

    print(f"Metadata table saved to {output_file}")

if __name__ == '__main__':
    main()
