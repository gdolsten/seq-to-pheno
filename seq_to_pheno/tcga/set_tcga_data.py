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
import sys
import requests
from Bio.Seq import Seq


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
            # Remove 'chr' prefix if present
            if chrom.startswith('chr'):
                chrom = chrom[3:]
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

def annotate_vcf(input_vcf: str, output_vcf: str, snpeff_jar: str, genome_version: str = 'GRCh37.75') -> None:
    """
    Annotates a VCF file using SnpEff.

    Parameters:
        input_vcf (str): Path to the input VCF file.
        output_vcf (str): Path to the output annotated VCF file.
        snpeff_jar (str): Path to the snpEff.jar file.
        genome_version (str): Genome version to use for annotation.

    Returns:
        None
    """
    try:
        command = [
            'java', '-Xmx4g', '-jar', snpeff_jar,
            '-v', genome_version,
            '-noLog',
            input_vcf
        ]
        logging.info(f"Annotating VCF file {input_vcf}")
        with open(output_vcf, 'w') as outfile:
            subprocess.run(command, stdout=outfile, stderr=subprocess.PIPE, check=True)
        logging.info(f"Annotated VCF saved to {output_vcf}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error annotating VCF file {input_vcf}: {e.stderr.decode().strip()}")

def count_variants_in_transcripts(vcf_file: str, interval_trees: dict) -> dict:
    """
    Counts the number of impactful variants overlapping each transcript in a VCF file.
    """
    variant_counts = defaultdict(int)
    try:
        # Suppress warnings from pysam/htslib
        save_stderr = sys.stderr
        sys.stderr = open(os.devnull, 'w')

        vcf = pysam.VariantFile(vcf_file)

        # Restore stderr
        sys.stderr.close()
        sys.stderr = save_stderr
    except Exception as e:
        logging.error(f"Failed to open VCF file {vcf_file}: {e}")
        return variant_counts

    for record in vcf:
        chrom = record.chrom
        # Remove 'chr' prefix if present
        if chrom.startswith('chr'):
            chrom = chrom[3:]
        pos = record.pos
        if chrom not in interval_trees:
            continue
        overlaps = interval_trees[chrom][pos]
        if not overlaps:
            continue

        # Parse the ANN field
        ann_field = record.info.get('ANN')
        if ann_field:
            for ann in ann_field:
                ann_parts = ann.split('|')
                # ANN fields are as follows:
                # Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO

                # We need Feature_ID (which is the transcript ID) and the Impact
                impact = ann_parts[2]
                feature_id = ann_parts[6]

                # Check if the feature_id matches any of our transcripts
                if feature_id in [interval.data for interval in overlaps]:
                    # Optionally, filter by impact level
                    if impact in ['HIGH', 'MODERATE']:
                        variant_counts[feature_id] += 1
        else:
            # No annotation available, skip this variant
            continue
    return variant_counts

def process_sample(args):
    """
    Processes a single sample to count variants per transcript.

    Parameters:
        args (tuple): A tuple containing sample_id, wgs_aliquot_id, interval_trees, variant_dir, snpeff_jar

    Returns:
        pd.DataFrame or None: DataFrame with transcript IDs as index and variant counts for this sample, or None if no data.
    """
    sample_id, wgs_aliquot_id, interval_trees, variant_dir, snpeff_jar = args
    variant_counts = defaultdict(int)

    # Paths to VCF files
    snv_vcf_pattern = os.path.join(variant_dir, 'snv_mnv', f"{wgs_aliquot_id}.consensus.*.somatic.snv_mnv.vcf.gz")
    indel_vcf_pattern = os.path.join(variant_dir, 'indel', f"{wgs_aliquot_id}.consensus.*.somatic.indel.vcf.gz")

    # Find actual VCF files using glob
    snv_vcf_files = glob.glob(snv_vcf_pattern)
    indel_vcf_files = glob.glob(indel_vcf_pattern)

    # Process SNV VCF
    if snv_vcf_files:
        snv_vcf = snv_vcf_files[0]
        annotated_snv_vcf = snv_vcf.replace('.vcf.gz', '.annotated.vcf')
        if not os.path.exists(annotated_snv_vcf):
            annotate_vcf(snv_vcf, annotated_snv_vcf, snpeff_jar)
        logging.info(f"Processing annotated SNV VCF for sample {sample_id}")
        snv_counts = count_variants_in_transcripts(annotated_snv_vcf, interval_trees)
        for k, v in snv_counts.items():
            variant_counts[k] += v
    else:
        logging.warning(f"SNV VCF not found for sample {sample_id}")

    # Process Indel VCF
    if indel_vcf_files:
        indel_vcf = indel_vcf_files[0]
        annotated_indel_vcf = indel_vcf.replace('.vcf.gz', '.annotated.vcf')
        if not os.path.exists(annotated_indel_vcf):
            annotate_vcf(indel_vcf, annotated_indel_vcf, snpeff_jar)
        logging.info(f"Processing annotated Indel VCF for sample {sample_id}")
        indel_counts = count_variants_in_transcripts(annotated_indel_vcf, interval_trees)
        for k, v in indel_counts.items():
            variant_counts[k] += v
    else:
        logging.warning(f"Indel VCF not found for sample {sample_id}")

    if variant_counts:
        # Create DataFrame for this sample
        df_sample = pd.DataFrame.from_dict(variant_counts, orient='index', columns=[sample_id])
        return df_sample
    else:
        logging.warning(f"No variants found for sample {sample_id}")
        return None

def process_all_samples(metadata_df: pd.DataFrame, interval_trees: dict, variant_dir: str, snpeff_jar: str) -> pd.DataFrame:
    """
    Processes all samples to count variants per transcript.

    Parameters:
        metadata_df (pd.DataFrame): Metadata DataFrame with samples to process.
        interval_trees (dict): Interval trees of transcripts.
        variant_dir (str): Directory containing variant files.
        snpeff_jar (str): Path to snpEff.jar file.

    Returns:
        pd.DataFrame: DataFrame with transcripts as rows and samples as columns.
    """
    data_frames = []
    for idx, row in metadata_df.iterrows():
        sample_id = row['aliquot_id']
        wgs_aliquot_id = row['wgs_aliquot_id']
        logging.info(f"Processing sample {sample_id}")
        args = (sample_id, wgs_aliquot_id, interval_trees, variant_dir, snpeff_jar)
        df_sample = process_sample(args)
        if df_sample is not None:
            data_frames.append(df_sample)

    # Concatenate all sample DataFrames
    if data_frames:
        variant_counts_df = pd.concat(data_frames, axis=1)
        # Fill missing values with 0
        variant_counts_df = variant_counts_df.fillna(0).astype(int)
        variant_counts_df.index.name = 'transcript_id'
        variant_counts_df = variant_counts_df.sort_index()
    else:
        variant_counts_df = pd.DataFrame()

    return variant_counts_df

def main():

    # Define constants or load them from a config file
    ENDPOINT_URL = 'https://object.genomeinformatics.org'
    DATA_DIR = 'seq_to_pheno/tcga/data'
    VARIANT_DIR = os.path.join(DATA_DIR, 'variants')
    TRANSCRIPT_DATA_S3 = 's3://icgc25k-open/PCAWG/transcriptome/transcript_expression/pcawg.rnaseq.transcript.expr.tpm.tsv.gz'
    METADATA_S3 = 's3://icgc25k-open/PCAWG/transcriptome/metadata/rnaseq.extended.metadata.aliquot_id.V4.tsv.gz'
    SNV_INDEL_S3 = 's3://icgc25k-open/PCAWG/consensus_snv_indel/final_consensus_snv_indel_passonly_icgc.public.tgz'
    SNPEFF_JAR = '/Users/harrison.reed/snpEff/snpEff.jar'  # Update with the correct path to your snpEff.jar
    GENOME_VERSION = 'GRCh37.75'      # Ensure this is the genome version you downloaded for snpEff

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
    metadata_df = metadata_df[metadata_df.has_consensus_data == True]
    metadata_df.to_csv(updated_metadata_file, sep='\t', index=False)
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
    variant_counts_df = process_all_samples(metadata_df, interval_trees, VARIANT_DIR, SNPEFF_JAR)

    # Save variant counts table
    variant_counts_file = os.path.join(DATA_DIR, 'variant_counts_per_transcript.tsv')
    variant_counts_df.to_csv(variant_counts_file, sep='\t', header=True)
    logging.info(f"Variant counts per transcript saved to {variant_counts_file}")


if __name__ == '__main__':
    main()