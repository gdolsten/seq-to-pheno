import pysam
import requests
import re
import os
import time
import logging
import pickle
from Bio.Seq import Seq
from concurrent.futures import ProcessPoolExecutor

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

# Initialize a cache for CDS sequences
if os.path.exists('cds_cache.pkl'):
    with open('cds_cache.pkl', 'rb') as f:
        transcript_cds_cache = pickle.load(f)
else:
    transcript_cds_cache = {}

session = requests.Session()

def fetch_transcript_cds(transcript_id):
    if transcript_id in transcript_cds_cache:
        return transcript_cds_cache[transcript_id]
    server = "https://rest.ensembl.org"
    ext = f"/sequence/id/{transcript_id}?type=cds"
    headers = {"Content-Type": "text/plain"}

    retries = 5
    backoff_factor = 1
    for attempt in range(retries):
        try:
            response = session.get(server + ext, headers=headers)
            if response.status_code == 429:
                retry_after = int(response.headers.get("Retry-After", backoff_factor))
                logging.warning(f"Rate limit exceeded. Retrying after {retry_after} seconds.")
                time.sleep(retry_after)
                backoff_factor *= 2
                continue
            elif not response.ok:
                logging.error(f"Error fetching CDS for {transcript_id}: {response.text}")
                return None
            else:
                sequence = response.text
                transcript_cds_cache[transcript_id] = sequence
                return sequence
        except requests.exceptions.RequestException as e:
            logging.error(f"RequestException fetching CDS for {transcript_id}: {e}")
            time.sleep(backoff_factor)
            backoff_factor *= 2
            continue
    logging.error(f"Failed to fetch CDS for {transcript_id} after {retries} attempts.")
    return None

def apply_variant_to_cds(sequence, hgvs_c):
    """
    Applies a variant to the CDS sequence based on the HGVS cDNA notation.

    Parameters:
        sequence (str): The reference CDS sequence.
        hgvs_c (str): The HGVS cDNA notation (e.g., 'c.123A>T').

    Returns:
        str: The mutated CDS sequence.
    """

    # Substitution (SNV)
    snv_pattern = r'c\.(\d+)([A-Z])>([A-Z])'

    # Deletion of a single nucleotide
    del_single_pattern = r'c\.(\d+)del([A-Z])?'

    # Deletion of a range
    del_range_pattern = r'c\.(\d+)_(\d+)del([A-Z]+)?'

    # Insertion
    ins_pattern = r'c\.(\d+)_(\d+)ins([A-Z]+)'

    # Deletion-Insertion (Indel)
    delins_single_pattern = r'c\.(\d+)delins([A-Z]+)'
    delins_range_pattern = r'c\.(\d+)_(\d+)delins([A-Z]+)'

    # Duplication
    dup_pattern = r'c\.(\d+)_(\d+)dup([A-Z]+)?'
    dup_single_pattern = r'c\.(\d+)dup([A-Z]+)?'

    # Inversion
    inv_pattern = r'c\.(\d+)_(\d+)inv'

    # Handle SNVs
    match = re.match(snv_pattern, hgvs_c)
    if match:
        pos = int(match.group(1)) - 1  # Convert to 0-based index
        ref_base = match.group(2)
        alt_base = match.group(3)
        # Validate reference base
        if sequence[pos] != ref_base:
            raise ValueError(f"Reference base mismatch at position {pos+1}: expected {ref_base}, found {sequence[pos]}")
        # Apply mutation
        mutated_sequence = sequence[:pos] + alt_base + sequence[pos+1:]
        return mutated_sequence

    # Handle single nucleotide deletions
    match = re.match(del_single_pattern, hgvs_c)
    if match:
        pos = int(match.group(1)) - 1
        del_base = match.group(2)
        if del_base:
            # Validate reference base
            if sequence[pos] != del_base:
                raise ValueError(f"Reference base mismatch at position {pos+1}: expected {del_base}, found {sequence[pos]}")
        # Apply deletion
        mutated_sequence = sequence[:pos] + sequence[pos+1:]
        return mutated_sequence

    # Handle range deletions
    match = re.match(del_range_pattern, hgvs_c)
    if match:
        start_pos = int(match.group(1)) - 1
        end_pos = int(match.group(2))
        del_bases = match.group(3)
        if del_bases:
            ref_seq = sequence[start_pos:end_pos]
            if ref_seq != del_bases:
                raise ValueError(f"Reference bases mismatch for deletion at positions {start_pos+1}-{end_pos}: expected {del_bases}, found {ref_seq}")
        # Apply deletion
        mutated_sequence = sequence[:start_pos] + sequence[end_pos:]
        return mutated_sequence

    # Handle insertions
    match = re.match(ins_pattern, hgvs_c)
    if match:
        pos1 = int(match.group(1)) - 1
        pos2 = int(match.group(2))
        ins_bases = match.group(3)
        # Apply insertion
        mutated_sequence = sequence[:pos2] + ins_bases + sequence[pos2:]
        return mutated_sequence

    # Handle single position delins
    match = re.match(delins_single_pattern, hgvs_c)
    if match:
        pos = int(match.group(1)) - 1
        ins_bases = match.group(2)
        # Delete one base and insert new bases
        mutated_sequence = sequence[:pos] + ins_bases + sequence[pos+1:]
        return mutated_sequence

    # Handle range delins
    match = re.match(delins_range_pattern, hgvs_c)
    if match:
        start_pos = int(match.group(1)) - 1
        end_pos = int(match.group(2))
        ins_bases = match.group(3)
        # Apply delins
        mutated_sequence = sequence[:start_pos] + ins_bases + sequence[end_pos:]
        return mutated_sequence

    # Handle single nucleotide duplications
    match = re.match(dup_single_pattern, hgvs_c)
    if match:
        pos = int(match.group(1)) - 1
        dup_base = match.group(2)
        if dup_base:
            # Validate reference base
            if sequence[pos] != dup_base:
                raise ValueError(f"Reference base mismatch at position {pos+1}: expected {dup_base}, found {sequence[pos]}")
            duplicated_sequence = dup_base
        else:
            duplicated_sequence = sequence[pos]
        # Apply duplication
        mutated_sequence = sequence[:pos+1] + duplicated_sequence + sequence[pos+1:]
        return mutated_sequence

    # Handle range duplications
    match = re.match(dup_pattern, hgvs_c)
    if match:
        start_pos = int(match.group(1)) - 1
        end_pos = int(match.group(2))
        dup_bases = match.group(3)
        if dup_bases:
            ref_seq = sequence[start_pos:end_pos]
            if ref_seq != dup_bases:
                raise ValueError(f"Reference bases mismatch for duplication at positions {start_pos+1}-{end_pos}: expected {dup_bases}, found {ref_seq}")
            duplicated_sequence = dup_bases
        else:
            duplicated_sequence = sequence[start_pos:end_pos]
        # Apply duplication
        mutated_sequence = sequence[:end_pos] + duplicated_sequence + sequence[end_pos:]
        return mutated_sequence

    # Handle inversions
    match = re.match(inv_pattern, hgvs_c)
    if match:
        start_pos = int(match.group(1)) - 1
        end_pos = int(match.group(2))
        ref_seq = sequence[start_pos:end_pos]
        inverted_sequence = ref_seq[::-1]  # Reverse the sequence
        # Apply inversion
        mutated_sequence = sequence[:start_pos] + inverted_sequence + sequence[end_pos:]
        return mutated_sequence

    # If mutation type not recognized
    raise NotImplementedError(f"Mutation type in HGVS notation '{hgvs_c}' not supported.")


def translate_cds_to_protein(cds_sequence):
    coding_dna = Seq(cds_sequence)
    protein_sequence = coding_dna.translate(to_stop=False)

    # Check for premature stop codons
    if '*' in protein_sequence[:-1]:  # Exclude the last character if it's a stop codon
        logging.warning("Premature stop codon detected in protein sequence.")

    return str(protein_sequence)

def save_protein_sequence_to_fasta(sequence, protein_id, output_file):
    with open(output_file, 'w') as f:
        f.write(f">{protein_id}\n")
        # Split sequence into lines of 60 characters
        for i in range(0, len(sequence), 60):
            f.write(sequence[i:i+60] + '\n')

def process_sample(args):
    vcf_file, output_dir, sample_id, wildtype_output_dir = args
    process_vcf(vcf_file, output_dir, sample_id, wildtype_output_dir)

def process_vcf_directory(vcf_dir, output_dir, wildtype_output_dir):
    tasks = []
    for filename in os.listdir(vcf_dir):
        if filename.endswith('.vcf') or filename.endswith('.vcf.gz'):
            sample_id = filename.replace('.combined.vcf', '').replace('.vcf', '').replace('.gz', '')
            vcf_file = os.path.join(vcf_dir, filename)
            tasks.append((vcf_file, output_dir, sample_id, wildtype_output_dir))
    return tasks

def process_vcf(vcf_file, output_dir, sample_id, wildtype_output_dir):
    vcf_in = pysam.VariantFile(vcf_file)
    transcript_variants = {}
    failed_variants = []

    for record in vcf_in:
        ann_field = record.info.get('ANN')
        if ann_field:
            for ann in ann_field:
                ann_parts = ann.split('|')
                if len(ann_parts) < 15:
                    logging.warning(f"Skipping malformed ANN entry: {ann}")
                    continue
                allele = ann_parts[0]
                annotation = ann_parts[1]
                impact = ann_parts[2]
                gene_name = ann_parts[3]
                gene_id = ann_parts[4]
                feature_type = ann_parts[5]
                transcript_id = ann_parts[6]
                transcript_biotype = ann_parts[7]
                hgvs_c = ann_parts[9]

                # Skip if HGVS.c is not available
                if not hgvs_c or hgvs_c == '':
                    continue

                # Proceed only if feature_type is 'transcript' and transcript_id starts with 'ENST'
                if feature_type != 'transcript' or not transcript_id.startswith('ENST'):
                    continue

                # Remove transcript version if present
                transcript_id = transcript_id.split('.')[0]

                # Filter for protein-coding transcripts
                if transcript_biotype != 'protein_coding':
                    continue

                # Filter by Annotation Impact
                acceptable_impacts = ['HIGH', 'MODERATE']
                if impact not in acceptable_impacts:
                    continue

                # Filter by Annotation
                acceptable_annotations = ['missense_variant', 'nonsense_variant', 'frameshift_variant',
                                          'splice_acceptor_variant', 'splice_donor_variant', 'start_lost', 'stop_lost', 'stop_gained']
                if annotation not in acceptable_annotations:
                    continue

                # Collect variants per transcript
                if transcript_id not in transcript_variants:
                    transcript_variants[transcript_id] = []
                transcript_variants[transcript_id].append((hgvs_c, record.pos))

    for transcript_id, variants in transcript_variants.items():
        # Fetch the reference CDS sequence
        cds_sequence = fetch_transcript_cds(transcript_id)
        if not cds_sequence:
            logging.warning(f"CDS sequence not found for transcript {transcript_id}")
            continue

        # Save wildtype protein sequence if not already saved
        wildtype_protein_file = os.path.join(wildtype_output_dir, f"{transcript_id}.fasta")
        if not os.path.exists(wildtype_protein_file):
            wildtype_protein_sequence = translate_cds_to_protein(cds_sequence)
            save_protein_sequence_to_fasta(wildtype_protein_sequence, transcript_id, wildtype_protein_file)
            logging.info(f"Saved wildtype protein sequence for {transcript_id} to {wildtype_protein_file}")

        protein_id = f"{sample_id}_{transcript_id}_mutated"
        output_file = os.path.join(output_dir, f"{protein_id}.fasta")
        
        if not os.path.exists(output_file):
          # Sort variants by cDNA position (descending)
          def extract_cdna_position(hgvs_c):
              match = re.match(r'c\.([\d+]+)', hgvs_c)
              if match:
                  return int(match.group(1))
              else:
                  logging.warning(f"Could not extract cDNA position from {hgvs_c}")
                  return float('inf')

          variants.sort(key=lambda x: extract_cdna_position(x[0]), reverse=True)

          mutated_cds_sequence = cds_sequence
          for hgvs_c, _ in variants:
              try:
                  mutated_cds_sequence = apply_variant_to_cds(mutated_cds_sequence, hgvs_c)
              except NotImplementedError as e:
                  logging.warning(f"Variant not implemented {hgvs_c} for transcript {transcript_id}: {e}")
                  failed_variants.append((transcript_id, hgvs_c, str(e)))
                  continue
              except Exception as e:
                  logging.error(f"Error applying variant {hgvs_c} to transcript {transcript_id}: {e}")
                  failed_variants.append((transcript_id, hgvs_c, str(e)))
                  continue
        
          # Translate the mutated CDS into protein sequence
          protein_sequence = translate_cds_to_protein(mutated_cds_sequence)

          # Save the mutated protein sequence to a FASTA file
          save_protein_sequence_to_fasta(protein_sequence, protein_id, output_file)
          logging.info(f"Saved mutated protein sequence for {transcript_id} to {output_file}")

    # Optionally, write failed variants to a log file
    if failed_variants:
        failed_variants_file = os.path.join(output_dir, f'{sample_id}_failed_variants.log')
        with open(failed_variants_file, 'w') as f:
            for transcript_id, hgvs_c, error_msg in failed_variants:
                f.write(f"{transcript_id}\t{hgvs_c}\t{error_msg}\n")
        logging.info(f"Failed variants logged to {failed_variants_file}")

    # Save the updated CDS cache
    with open('cds_cache.pkl', 'wb') as f:
        pickle.dump(transcript_cds_cache, f)

if __name__ == '__main__':
    # Load cache from file if exists
    if os.path.exists('cds_cache.pkl'):
        with open('cds_cache.pkl', 'rb') as f:
            transcript_cds_cache = pickle.load(f)
    else:
        transcript_cds_cache = {}

    vcf_dir = 'seq_to_pheno/tcga/data/variants/vcf'
    output_dir = 'seq_to_pheno/tcga/data/variants/mutated_proteins'
    wildtype_output_dir = 'seq_to_pheno/tcga/data/variants/wildtype_proteins'
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(wildtype_output_dir, exist_ok=True)

    tasks = process_vcf_directory(vcf_dir, output_dir, wildtype_output_dir)

    # Start the multiprocessing executor
    with ProcessPoolExecutor(max_workers=4) as executor:
        executor.map(process_sample, tasks)

    # Save the updated CDS cache
    with open('cds_cache.pkl', 'wb') as f:
        pickle.dump(transcript_cds_cache, f)
