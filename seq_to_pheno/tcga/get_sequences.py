import pysam
import requests
import re
import os
import time
import logging
from Bio.Seq import Seq

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

# Initialize a cache for CDS sequences
transcript_cds_cache = {}


def fetch_transcript_cds(transcript_id):
    if transcript_id in transcript_cds_cache:
        return transcript_cds_cache[transcript_id]
    server = "https://rest.ensembl.org"
    ext = f"/sequence/id/{transcript_id}?type=cds"
    headers = {"Content-Type": "text/plain"}

    while True:
        response = requests.get(server + ext, headers=headers)
        if response.status_code == 429:
            # Too Many Requests, wait and retry
            retry_after = int(response.headers.get("Retry-After", 1))
            logging.warning(f"Rate limit exceeded. Retrying after {retry_after} seconds.")
            time.sleep(retry_after)
            continue
        elif not response.ok:
            logging.error(f"Error fetching CDS for {transcript_id}: {response.text}")
            return None
        else:
            sequence = response.text
            transcript_cds_cache[transcript_id] = sequence
            return sequence


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


def process_vcf(vcf_file, output_dir):
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

        # Sort variants by cDNA position (descending)
        try:
            variants.sort(key=lambda x: int(re.match(r'c\.(\d+)', x[0]).group(1)), reverse=True)
        except Exception as e:
            logging.error(f"Error sorting variants for transcript {transcript_id}: {e}")
            continue

        mutated_cds_sequence = cds_sequence
        for hgvs_c, _ in variants:
            try:
                mutated_cds_sequence = apply_variant_to_cds(mutated_cds_sequence, hgvs_c)
            except Exception as e:
                logging.error(f"Error applying variant {hgvs_c} to transcript {transcript_id}: {e}")
                failed_variants.append((transcript_id, hgvs_c, str(e)))
                continue

        # Translate the mutated CDS into protein sequence
        protein_sequence = translate_cds_to_protein(mutated_cds_sequence)

        # Save the mutated protein sequence to a FASTA file
        protein_id = f"{transcript_id}_mutated"
        output_file = os.path.join(output_dir, f"{protein_id}.fasta")
        save_protein_sequence_to_fasta(protein_sequence, protein_id, output_file)
        logging.info(f"Saved mutated protein sequence for {transcript_id} to {output_file}")

    # Optionally, write failed variants to a log file
    if failed_variants:
        failed_variants_file = os.path.join(output_dir, 'failed_variants.log')
        with open(failed_variants_file, 'w') as f:
            for transcript_id, hgvs_c, error_msg in failed_variants:
                f.write(f"{transcript_id}\t{hgvs_c}\t{error_msg}\n")
        logging.info(f"Failed variants logged to {failed_variants_file}")


# Example usage
vcf_file = 'tcga/data/variants/snv_mnv/d8c2b4b2-e12b-43d2-bafc-87b29f027797.consensus.20160830.somatic.snv_mnv.annotated.vcf'
output_dir = 'tcga/data/variants/snv_mnv/mutated_transcripts'
os.makedirs(output_dir, exist_ok=True)

process_vcf(vcf_file, output_dir)
