import argparse
import os
from seq_to_pheno.longevity.sequence_object.SequenceObject import SpeciesSequenceObject
from seq_to_pheno.hug.utils import DatasetPusher, DatasetCard
from typing import Dict, Any

def prepare_card_data(output_file: str, filtered_mapped_ortholog_df: Any) -> Dict[str, Any]:
    dataset_size = os.path.getsize(output_file)
    num_examples = len(filtered_mapped_ortholog_df)
    
    card_data = {
        "license": "CC-BY-4.0",
        "size_category": "100M<n<1G" if dataset_size < 1e9 else "1G<n<10G",
        "task_category": "other",
        "pretty_name": "ZoonomiaOrthologs",
        "num_bytes": dataset_size,
        "num_examples": num_examples,
        "download_size": dataset_size,
        "dataset_size": dataset_size,
        "data_file_path": os.path.basename(output_file)
    }

    return card_data

def main():
    # Argument parser setup
    parser = argparse.ArgumentParser(description="Extract protein sequences and optionally publish to Hugging Face Hub.")
    
    # Required arguments for sequence extraction
    parser.add_argument('--input_folder', type=str, required=True, help='Path to the folder containing the protein sequence data')
    parser.add_argument('--input_file', type=str, required=True, help='Protein sequence file (e.g., protein_sequence_df.tsv)')
    parser.add_argument('--output_folder', type=str, required=True, help='Folder to save the extracted protein sequences')
    parser.add_argument('--output_file', type=str, required=True, help='Output FASTA file for protein sequences')
    parser.add_argument('--gene', type=str, default=None, help='Gene symbol to extract sequences for (default: extract all genes)')

    # Optional arguments for Hugging Face Hub upload
    parser.add_argument('--publish', action='store_true', help='Whether to publish the dataset to Hugging Face Hub')
    parser.add_argument('--repo_name', type=str, help='Name of the dataset repository on Hugging Face Hub')
    parser.add_argument('--hf_token', type=str, default=None, help='Hugging Face token for authentication (or set via HF_TOKEN env variable)')
    parser.add_argument('--card_template', type=str, default='seq_to_pheno/hug/zoonomia_fasta_repo_template/README.md', help='Path to the dataset card template')

    args = parser.parse_args()

    # Ensure output directory exists
    os.makedirs(args.output_folder, exist_ok=True)

    # Load species sequence data
    species = SpeciesSequenceObject.load(args.input_folder, args.input_file)

    # Determine output file path
    output_path = os.path.join(args.output_folder, args.output_file)

    # Extract sequences
    if args.gene:
        # Extract sequences for a specific gene
        species.to_fasta(args.output_folder, output_path, args.gene)
        print(f"Protein sequences for gene {args.gene} extracted to {output_path}.")
    else:
        # Extract sequences for all genes
        species.to_fasta(args.output_folder, output_path)
        print(f"All protein sequences extracted to {output_path}.")

    # If publish flag is set, push dataset to Hugging Face Hub
    if args.publish:
        if not args.repo_name:
            raise ValueError("You must provide a repository name (--repo_name) when publishing the dataset to Hugging Face Hub.")

        # Set up the dataset card and token
        token = args.hf_token or os.environ.get("HF_TOKEN")
        if not token:
            raise ValueError("Hugging Face token is required. Set it via --hf_token or HF_TOKEN environment variable.")

        # Prepare dataset card
        card_data = prepare_card_data(output_file=output_path, filtered_mapped_ortholog_df=species.sequence_df)
        dataset_card = DatasetCard(args.card_template)

        # Create and push dataset to the Hugging Face Hub
        pusher = DatasetPusher(folder=args.output_folder, name=args.repo_name, token=token)
        pusher.push_to_hub(dataset_card, card_data)

        print(f"Dataset successfully published to Hugging Face Hub under the repository: {args.repo_name}")

if __name__ == "__main__":
    main()
