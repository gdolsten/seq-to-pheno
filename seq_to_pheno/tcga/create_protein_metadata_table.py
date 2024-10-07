import os
import pandas as pd
import argparse

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

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Create protein sequences metadata table.')
    parser.add_argument('--include-sequences', action='store_true',
                        help='Include actual protein sequences in the output instead of file paths.')
    args = parser.parse_args()

    # Paths to directories
    mutated_proteins_dir = 'seq_to_pheno/tcga/data/variants/mutated_proteins'
    wildtype_proteins_dir = 'seq_to_pheno/tcga/data/variants/wildtype_proteins'
    metadata_file = 'seq_to_pheno/tcga/data/rnaseq.extended.metadata.aliquot_id.V4.tsv.gz'
    output_file = 'seq_to_pheno/tcga/data/protein_sequences_metadata.tsv'

    # Read the metadata file into a DataFrame
    metadata_df = pd.read_csv(metadata_file, sep='\t', compression='gzip', low_memory=False)

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
        # Ensure to handle missing data
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
            row['mutated_protein_sequence'] = mutated_sequence
            row['wildtype_protein_sequence'] = wildtype_sequence
        else:
            # Include the file paths
            row['mutated_protein_path'] = mutated_protein_path
            row['wildtype_protein_path'] = wildtype_protein_path

        for col_name, meta_col in columns_to_extract.items():
            value = sample_metadata.get(meta_col, '')
            row[col_name] = value

        data_rows.append(row)

    # Create a DataFrame from the collected rows
    result_df = pd.DataFrame(data_rows)

    # Save the DataFrame to a TSV file
    result_df.to_csv(output_file, index=False, sep='\t')

    print(f"Metadata table saved to {output_file}")

if __name__ == '__main__':
    main()
