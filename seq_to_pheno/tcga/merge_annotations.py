import os
import subprocess
import logging
import re

# Compress and index VCFs if necessary
def compress_and_index_vcf(vcf_file):
    if not vcf_file.endswith('.vcf.gz'):
        # Compress the VCF
        compressed_vcf = vcf_file + '.gz'
        logging.info(f'Compressing {vcf_file} to {compressed_vcf}')
        with open(compressed_vcf, 'wb') as f_out:
            subprocess.run(['bgzip', '-@', '4', '-c', vcf_file], stdout=f_out, check=True)
        # Index the compressed VCF
        logging.info(f'Indexing {compressed_vcf}')
        subprocess.run(['tabix', '-p', 'vcf', compressed_vcf], check=True)
        return compressed_vcf
    else:
        # Ensure the index exists
        if not os.path.exists(vcf_file + '.tbi'):
            logging.info(f'Indexing {vcf_file}')
            subprocess.run(['tabix', '-p', 'vcf', vcf_file], check=True)
        return vcf_file


# Configure logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

# Directories containing the annotated VCF files
indel_dir = 'seq_to_pheno/tcga/data/variants/indel/'
snv_mnv_dir = 'seq_to_pheno/tcga/data/variants/snv_mnv/'

# Output directory
output_dir = 'seq_to_pheno/tcga/data/variants/vcf/'
os.makedirs(output_dir, exist_ok=True)

# sample-name regex
pattern = '(.*).consensus.*'

# Get a list of sample names
sample_names = set()
indel_names = {}
snv_mnv_names = {}
for filename in os.listdir(indel_dir):
    if filename.endswith('.annotated.vcf') or filename.endswith('.annotated.vcf.gz'):
        sample_name = re.match(pattern, filename).group(1)
        sample_names.add(sample_name)
        indel_names[sample_name] = filename

for filename in os.listdir(snv_mnv_dir):
    if filename.endswith('.annotated.vcf') or filename.endswith('.annotated.vcf.gz'):
        sample_name = re.match(pattern, filename).group(1)
        snv_mnv_names[sample_name] = filename

# Process each sample
for sample_name in sample_names:
    indel_vcf = compress_and_index_vcf(os.path.join(indel_dir, indel_names[sample_name]))
    snv_vcf = compress_and_index_vcf(os.path.join(snv_mnv_dir, snv_mnv_names.get(sample_name)))
    combined_vcf = os.path.join(output_dir, f'{sample_name}.combined.vcf')

    # Check if both VCF files exist
    if not os.path.exists(indel_vcf):
        logging.warning(f'Indel VCF not found for sample {sample_name}: {indel_vcf}')
        continue

    if not snv_vcf or not os.path.exists(snv_vcf):
        logging.warning(f'SNV/MNV VCF not found for sample {sample_name}: {snv_vcf}')
        continue

    # Decide output type based on input file extensions
    if indel_vcf.endswith('.vcf.gz') or snv_vcf.endswith('.vcf.gz'):
        # Use compressed output
        combined_vcf += '.gz'
        output_type = 'z'
    else:
        output_type = 'v'

    # Use bcftools merge to combine the VCF files
    cmd = [
        'bcftools', 'merge',
        '--merge', 'all',
        '--output', combined_vcf,
        '--output-type', output_type,
        '--threads', '8',
        indel_vcf, snv_vcf
    ]

    logging.info(f'Combining VCFs for sample {sample_name}')
    try:
        subprocess.run(cmd, check=True)
        logging.info(f'Combined VCF written to {combined_vcf}')
        # Index the combined VCF if compressed
        if output_type == 'z':
            subprocess.run(['bcftools', 'index', combined_vcf], check=True)
            logging.info(f'Indexed combined VCF: {combined_vcf}')
    except subprocess.CalledProcessError as e:
        logging.error(f'Error combining VCFs for sample {sample_name}: {e}')
