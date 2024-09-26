import os
import subprocess
import logging
import pandas as pd
from typing import Set
import gzip
from collections import defaultdict
import pysam
from intervaltree import IntervalTree
import glob


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


def parse_gtf(gtf_file: str) -> pd.DataFrame:
    """
    Parses a GTF file to extract transcript genomic coordinates.

    Parameters:
        gtf_file (str): Path to the GTF file.

    Returns:
        pd.DataFrame: DataFrame with columns ['transcript_id', 'chr', 'start', 'end', 'strand'].
    """
    transcripts = []
    with gzip.open(gtf_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if fields[2] != 'transcript':
                continue
            chrom = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            info_field = fields[8]
            attributes = {item.strip().split(' ')[0]: item.strip().split(' ')[1].replace('"', '') for item in info_field.strip().split(';') if item}
            transcript_id = attributes.get('transcript_id')
            if transcript_id:
                transcripts.append({
                    'transcript_id': transcript_id,
                    'chr': chrom,
                    'start': start,
                    'end': end,
                    'strand': strand
                })
    return pd.DataFrame(transcripts)


def build_interval_tree(transcript_df: pd.DataFrame) -> dict:
    """
    Builds an interval tree for each chromosome from the transcript DataFrame.

    Parameters:
        transcript_df (pd.DataFrame): DataFrame with transcript coordinates.

    Returns:
        dict: A dictionary with chromosome names as keys and IntervalTrees as values.
    """
    interval_trees = defaultdict(IntervalTree)
    for _, row in transcript_df.iterrows():
        chrom = row['chr']
        interval_trees[chrom].addi(row['start'], row['end'], row['transcript_id'])
    return interval_trees

def count_variants_in_transcripts(vcf_file: str, interval_trees: dict) -> dict:
    """
    Counts the number of variants overlapping each transcript in a VCF file.

    Parameters:
        vcf_file (str): Path to the VCF file.
        interval_trees (dict): Interval trees for each chromosome.

    Returns:
        dict: A dictionary with transcript IDs as keys and variant counts as values.
    """
    variant_counts = defaultdict(int)
    vcf = pysam.VariantFile(vcf_file)
    for record in vcf:
        chrom = record.chrom
        pos = record.pos
        if chrom not in interval_trees:
            continue
        overlaps = interval_trees[chrom][pos]
        for overlap in overlaps:
            transcript_id = overlap.data
            variant_counts[transcript_id] += 1
    return variant_counts

def process_all_samples(metadata_df: pd.DataFrame, interval_trees: dict, variant_dir: str) -> pd.DataFrame:
    """
    Processes all samples to count variants per transcript.

    Parameters:
        metadata_df (pd.DataFrame): Metadata DataFrame with samples to process.
        interval_trees (dict): Interval trees of transcripts.
        variant_dir (str): Directory containing variant files.

    Returns:
        pd.DataFrame: DataFrame with transcripts as rows and samples as columns.
    """
    data_frames = []
    
    for idx, row in metadata_df.iterrows():
        sample_id = row['aliquot_id']
        wgs_aliquot_id = row['wgs_aliquot_id']
        logging.info(f"Processing sample {sample_id}")
        
        # Paths to VCF files
        snv_vcf = os.path.join(variant_dir, 'snv_mnv', f"{wgs_aliquot_id}.consensus.*.somatic.snv_mnv.vcf.gz")
        indel_vcf = os.path.join(variant_dir, 'indel', f"{wgs_aliquot_id}.consensus.*.somatic.indel.vcf.gz")
        
        variant_counts = defaultdict(int)
        
        # Process SNV VCF
        if glob.glob(snv_vcf):
            snv_vcf_file = glob.glob(snv_vcf)[0]
            logging.info(f"Processing SNV VCF for sample {sample_id}")
            snv_counts = count_variants_in_transcripts(snv_vcf_file, interval_trees)
            for k, v in snv_counts.items():
                variant_counts[k] += v
        else:
            logging.warning(f"SNV VCF not found for sample {sample_id}")
        
        # Process Indel VCF
        if glob.glob(indel_vcf):
            indel_vcf_file = glob.glob(indel_vcf)[0]
            logging.info(f"Processing Indel VCF for sample {sample_id}")
            indel_counts = count_variants_in_transcripts(indel_vcf_file, interval_trees)
            for k, v in indel_counts.items():
                variant_counts[k] += v
        else:
            logging.warning(f"Indel VCF not found for sample {sample_id}")
        
        # Create DataFrame for this sample
        df_sample = pd.DataFrame.from_dict(variant_counts, orient='index', columns=[sample_id])
        data_frames.append(df_sample)
    
    # Concatenate all sample DataFrames
    variant_counts_df = pd.concat(data_frames, axis=1)
    # Fill missing values with 0
    variant_counts_df = variant_counts_df.fillna(0).astype(int)
    variant_counts_df.index.name = 'transcript_id'
    variant_counts_df = variant_counts_df.sort_index()
    
    return variant_counts_df

def main():

    # Define constants or load them from a config file
    ENDPOINT_URL = 'https://object.genomeinformatics.org'
    DATA_DIR = 'seq_to_pheno/tcga/data'
    VARIANT_DIR = os.path.join(DATA_DIR, 'variants')
    TRANSCRIPT_DATA_S3 = 's3://icgc25k-open/PCAWG/transcriptome/transcript_expression/pcawg.rnaseq.transcript.expr.tpm.tsv.gz'
    METADATA_S3 = 's3://icgc25k-open/PCAWG/transcriptome/metadata/rnaseq.extended.metadata.aliquot_id.V4.tsv.gz'
    SNV_INDEL_S3 = 's3://icgc25k-open/PCAWG/consensus_snv_indel/final_consensus_snv_indel_passonly_icgc.public.tgz'

    # Create data directories if they don't exist
    os.makedirs(DATA_DIR, exist_ok=True)
    os.makedirs(VARIANT_DIR, exist_ok=True)

    # Download files
    # download_file(TRANSCRIPT_DATA_S3, DATA_DIR, ENDPOINT_URL)
    # download_file(METADATA_S3, DATA_DIR, ENDPOINT_URL)
    # download_file(SNV_INDEL_S3, DATA_DIR, ENDPOINT_URL)

    # Extract SNV and Indel data
    # extract_tar_gz(os.path.join(DATA_DIR, 'final_consensus_snv_indel_passonly_icgc.public.tgz'), VARIANT_DIR)

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
    
    # Load Ensembl GTF file
    gtf_file = os.path.join(DATA_DIR, 'Homo_sapiens.GRCh37.87.gtf.gz')
    if not os.path.exists(gtf_file):
        logging.info("Downloading Ensembl GTF file")
        gtf_url = 'ftp://ftp.ensembl.org/pub/grch37/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz'
        download_file(gtf_url, DATA_DIR, '')
    
    logging.info("Parsing GTF file to get transcript coordinates")
    transcript_df = parse_gtf(gtf_file)
    
    logging.info("Building interval trees for transcripts")
    interval_trees = build_interval_tree(transcript_df)
    
    # Process samples to get variant counts per transcript
    logging.info("Processing samples to count variants per transcript")
    variant_counts_df = process_all_samples(metadata_df, interval_trees, VARIANT_DIR)
    
    # Save variant counts table
    variant_counts_file = os.path.join(DATA_DIR, 'variant_counts_per_transcript.tsv')
    variant_counts_df.to_csv(variant_counts_file, sep='\t', header=True)
    logging.info(f"Variant counts per transcript saved to {variant_counts_file}")
    
if __name__ == '__main__':
    main()
