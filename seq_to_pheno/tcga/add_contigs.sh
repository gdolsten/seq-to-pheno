#!/bin/bash

# Directories
# VCF_DIR="seq_to_pheno/tcga/data/variants/old_indel"    # Replace with the path to your VCF files
# OUTPUT_DIR="seq_to_pheno/tcga/data/variants/indel"    # Replace with the desired output directory
# CONTIGS_FILE="contigs.txt"           # Path to your contigs.txt file

# mkdir -p "$OUTPUT_DIR"

for vcf_file in "$VCF_DIR"/*.vcf.gz; do
    vcf_basename=$(basename "$vcf_file")
    output_file="$OUTPUT_DIR/$vcf_basename"

    echo "Processing $vcf_basename"

    # Add contig lines to the VCF header
    bcftools annotate --header-lines "$CONTIGS_FILE" -O z -o "$output_file" "$vcf_file"

    # Index the updated VCF file
    tabix -p vcf "$output_file"
done

VCF_DIR="seq_to_pheno/tcga/data/variants/old_snv_mnv"    # Replace with the path to your VCF files
OUTPUT_DIR="seq_to_pheno/tcga/data/variants/snv_mnv"    # Replace with the desired output directory
CONTIGS_FILE="contigs.txt"           # Path to your contigs.txt file

mkdir -p "$OUTPUT_DIR"

for vcf_file in "$VCF_DIR"/*.vcf.gz; do
    vcf_basename=$(basename "$vcf_file")
    output_file="$OUTPUT_DIR/$vcf_basename"

    echo "Processing $vcf_basename"

    # Add contig lines to the VCF header
    bcftools annotate --header-lines "$CONTIGS_FILE" -O z -o "$output_file" "$vcf_file"

    # Index the updated VCF file
    tabix -p vcf "$output_file"
done