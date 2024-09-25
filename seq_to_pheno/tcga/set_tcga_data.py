import os
import subprocess
import logging
import pandas as pd
from typing import Set

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

def download_file(s3_path: str, local_dir: str, endpoint_url: str) -> None:
    """
    Downloads a file from an S3 bucket using AWS CLI.

    Parameters:
        s3_path (str): The S3 path to the file.
        local_dir (str): The local directory where the file will be saved.
        endpoint_url (str): The endpoint URL for the S3 service.

    Returns:
        None
    """
    try:
        command = [
            'aws', 's3', 'cp', s3_path, local_dir,
            '--endpoint-url', endpoint_url, '--no-sign-request'
        ]
        logging.info(f"Downloading {s3_path} to {local_dir}")
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Failed to download {s3_path}: {e}")

def extract_tar_gz(file_path: str, extract_path: str) -> None:
    """
    Extracts a .tar.gz file to the specified directory.

    Parameters:
        file_path (str): The path to the .tar.gz file.
        extract_path (str): The directory to extract files into.

    Returns:
        None
    """
    try:
        logging.info(f"Extracting {file_path} to {extract_path}")
        subprocess.run(['tar', '-zxvf', file_path, '-C', extract_path], check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Failed to extract {file_path}: {e}")

def get_wgs_aliquot_ids_from_filenames(directory: str) -> Set[str]:
    """
    Extracts wgs_aliquot_ids from filenames in the given directory.

    Parameters:
        directory (str): The directory containing the variant files.

    Returns:
        Set[str]: A set of wgs_aliquot_ids extracted from the filenames.
    """
    wgs_aliquot_ids = set()
    for filename in os.listdir(directory):
        if filename.endswith('.vcf.gz'):
            # The filename format is {wgs_aliquot_id}.consensus.*.vcf.gz
            wgs_aliquot_id = filename.split('.')[0]
            wgs_aliquot_ids.add(wgs_aliquot_id)
    logging.info(f"Found {len(wgs_aliquot_ids)} wgs_aliquot_ids in {directory}")
    return wgs_aliquot_ids

def add_has_consensus_data_column(metadata_df: pd.DataFrame, wgs_ids_set: Set[str]) -> pd.DataFrame:
    """
    Adds a 'has_consensus_data' column to the metadata DataFrame.

    Parameters:
        metadata_df (pd.DataFrame): The original metadata DataFrame.
        wgs_ids_set (Set[str]): Set of wgs_aliquot_ids that have consensus data.

    Returns:
        pd.DataFrame: The updated metadata DataFrame with the new column.
    """
    metadata_df['has_consensus_data'] = metadata_df['wgs_aliquot_id'].isin(wgs_ids_set)
    logging.info("Added 'has_consensus_data' column to metadata DataFrame")
    return metadata_df

def main():

    # Define constants or load them from a config file
    ENDPOINT_URL = 'https://object.genomeinformatics.org'
    DATA_DIR = 'data'
    VARIANT_DIR = os.path.join(DATA_DIR, 'variants')
    TRANSCRIPT_DATA_S3 = 's3://icgc25k-open/PCAWG/transcriptome/transcript_expression/pcawg.rnaseq.transcript.expr.tpm.tsv.gz'
    METADATA_S3 = 's3://icgc25k-open/PCAWG/transcriptome/metadata/rnaseq.extended.metadata.aliquot_id.V4.tsv.gz'
    SNV_INDEL_S3 = 's3://icgc25k-open/PCAWG/consensus_snv_indel/final_consensus_snv_indel_passonly_icgc.public.tgz'

    # Create data directories if they don't exist
    os.makedirs(DATA_DIR, exist_ok=True)
    os.makedirs(VARIANT_DIR, exist_ok=True)

    # Download files
    download_file(TRANSCRIPT_DATA_S3, DATA_DIR, ENDPOINT_URL)
    download_file(METADATA_S3, DATA_DIR, ENDPOINT_URL)
    download_file(SNV_INDEL_S3, DATA_DIR, ENDPOINT_URL)

    # Extract SNV and Indel data
    extract_tar_gz(os.path.join(DATA_DIR, 'final_consensus_snv_indel_passonly_icgc.public.tgz'), VARIANT_DIR)

    # Paths to variant directories
    indel_dir = os.path.join(VARIANT_DIR, 'indel')
    snv_mnv_dir = os.path.join(VARIANT_DIR, 'snv_mnv')

    # Get wgs_aliquot_ids from filenames in both directories
    indel_wgs_ids = get_wgs_aliquot_ids_from_filenames(indel_dir)
    snv_wgs_ids = get_wgs_aliquot_ids_from_filenames(snv_mnv_dir)

    # Union of wgs_aliquot_ids from both directories
    all_wgs_ids = indel_wgs_ids.union(snv_wgs_ids)
    logging.info(f"Total unique wgs_aliquot_ids with consensus data: {len(all_wgs_ids)}")

    # Load metadata
    metadata_file = os.path.join(DATA_DIR, 'rnaseq.extended.metadata.aliquot_id.V4.tsv.gz')
    logging.info(f"Loading metadata from {metadata_file}")
    metadata_df = pd.read_csv(metadata_file, sep='\t', compression='gzip', low_memory=False)

    # Add 'has_consensus_data' column
    metadata_df = add_has_consensus_data_column(metadata_df, all_wgs_ids)

    # Save the updated metadata to a new file
    updated_metadata_file = os.path.join(DATA_DIR, 'metadata_with_consensus_data.tsv')
    
    num_samples_with_data = metadata_df['has_consensus_data'].sum()
    logging.info(f"Number of samples with consensus data: {num_samples_with_data}")
    
    # Filter metadata to include only samples with consensus data
    metadata_df[metadata_df.has_consensus_data == True].to_csv(updated_metadata_file, sep='\t', index=False)
    logging.info(f"Updated metadata saved to {updated_metadata_file}")
    
if __name__ == '__main__':
    main()
