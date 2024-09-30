import pysam
import requests
import re
import os
from Bio.Seq import Seq


def fetch_transcript_cds(transcript_id):
    server = "https://rest.ensembl.org"
    ext = f"/sequence/id/{transcript_id}?type=cds"
    headers = {"Content-Type": "text/plain"}
    response = requests.get(server + ext, headers=headers)
    if not response.ok:
        response.raise_for_status()
        return None
    sequence = response.text
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

    # Regular expression patterns for different mutation types
    snv_pattern = r'c\.(\d+)([A-Z])>([A-Z])'
    del_pattern = r'c\.(\d+)del([A-Z]+)?'
    ins_pattern = r'c\.(\d+)_(\d+)ins([A-Z]+)'
    delins_pattern = r'c\.(\d+)delins([A-Z]+)'
    delins_range_pattern = r'c\.(\d+)_(\d+)delins([A-Z]+)'

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

    # Handle deletions
    match = re.match(del_pattern, hgvs_c)
    if match:
        pos = int(match.group(1)) - 1
        del_bases = match.group(2)
        if del_bases:
            end_pos = pos + len(del_bases)
            if sequence[pos:end_pos] != del_bases:
                raise ValueError(f"Reference bases mismatch for deletion at positions {pos+1}-{end_pos}")
        else:
            end_pos = pos + 1
        mutated_sequence = sequence[:pos] + sequence[end_pos:]
        return mutated_sequence

    # Handle insertions
    match = re.match(ins_pattern, hgvs_c)
    if match:
        pos1 = int(match.group(1)) - 1
        pos2 = int(match.group(2))
        ins_bases = match.group(3)
        mutated_sequence = sequence[:pos2] + ins_bases + sequence[pos2:]
        return mutated_sequence

    # Handle delins (single position)
    match = re.match(delins_pattern, hgvs_c)
    if match:
        pos = int(match.group(1)) - 1
        ins_bases = match.group(2)
        # Delete one base and insert new bases
        mutated_sequence = sequence[:pos] + ins_bases + sequence[pos+1:]
        return mutated_sequence

    # Handle delins (range)
    match = re.match(delins_range_pattern, hgvs_c)
    if match:
        pos1 = int(match.group(1)) - 1
        pos2 = int(match.group(2))
        ins_bases = match.group(3)
        mutated_sequence = sequence[:pos1] + ins_bases + sequence[pos2:]
        return mutated_sequence

    # If mutation type not recognized
    raise NotImplementedError(f"Mutation type in HGVS notation '{hgvs_c}' not supported.")


def translate_cds_to_protein(cds_sequence):
    """
    Translates a CDS sequence into a protein sequence.

    Parameters:
        cds_sequence (str): The CDS nucleotide sequence.

    Returns:
        str: The protein amino acid sequence.
    """

    # Create a Seq object
    coding_dna = Seq(cds_sequence)

    # Translate the sequence
    protein_sequence = coding_dna.translate(to_stop=False)

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

    for record in vcf_in:
        ann_field = record.info.get('ANN')
        if ann_field:
            for ann in ann_field:
                ann_parts = ann.split('|')
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
            continue

        # Sort variants by position (descending)
        variants.sort(key=lambda x: x[1], reverse=True)

        mutated_cds_sequence = cds_sequence
        for hgvs_c, _ in variants:
            try:
                mutated_cds_sequence = apply_variant_to_cds(mutated_cds_sequence, hgvs_c)
            except Exception as e:
                print(f"Error applying variant {hgvs_c} to transcript {transcript_id}: {e}")
                continue

        # Translate the mutated CDS into protein sequence
        protein_sequence = translate_cds_to_protein(mutated_cds_sequence)

        # Save the mutated protein sequence to a FASTA file
        protein_id = f"{transcript_id}_mutated"
        output_file = os.path.join(output_dir, f"{protein_id}.fasta")
        save_protein_sequence_to_fasta(protein_sequence, protein_id, output_file)
        print(f"Saved mutated protein sequence for {transcript_id} to {output_file}")
        
        
# Example usage
vcf_file = 'tcga/data/variants/snv_mnv/d8c2b4b2-e12b-43d2-bafc-87b29f027797.consensus.20160830.somatic.snv_mnv.annotated.vcf'
output_dir = 'tcga/data/variants/snv_mnv/mutated_transcripts'
os.makedirs(output_dir, exist_ok=True)

process_vcf(vcf_file, output_dir)