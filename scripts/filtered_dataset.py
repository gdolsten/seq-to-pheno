import os
import argparse
from seq_to_pheno.hug.utils import DatasetCard, DatasetPusher
from seq_to_pheno.longevity.data.zoonomia.download_sequences_zoonomia import  fetch_directory_listing , download_files, filter_long_protein_sequences, filter_multimap_protein_sequences, count_mapped_orthologs
from datetime import datetime

import os
import argparse
from seq_to_pheno.hug.utils import DatasetCard, DatasetPusher
from datetime import datetime

def main():
    # Argument parser setup
    parser = argparse.ArgumentParser(description="Process and optionally publish a filtered Zoonomia dataset.")
    parser.add_argument('--folder', type=str, required=True, help='Path to the dataset folder')
    parser.add_argument('--output_folder', type=str, default='filtered_zoonomia_dataset', help='Path to save filtered dataset')
    parser.add_argument('--template_path', type=str, required=True, help='Path to the dataset README template file')
    parser.add_argument('--token', type=str, help='Hugging Face token for authentication (optional if --publish is not set)')
    parser.add_argument('--max_length', type=int, default=1000, help='Max protein sequence length for filtering')
    parser.add_argument('--max_orthologs', type=int, default=20, help='Max number of orthologs per protein for filtering')
    parser.add_argument('--publish', action='store_true', help='Publish the dataset to Hugging Face Hub')
    parser.add_argument('--animal_families', nargs='+', default=[
        'Afrotheria', 'Carnivora', 'Chiroptera', 'Dermoptera', 'Eulipotyphla', 'Lagomorpha', 'Metatheria', 
        'Perissodactyla', 'Pholidota', 'Primates', 'Prototheria', 'Rodentia', 'Ruminantia', 'Scandentia', 
        'Suina', 'Tylopoda', 'Whippomorpha', 'Xenarthra'
    ], help='List of animal families to include in the dataset')
    
    args = parser.parse_args()

    BASE_URL = "http://genome.senckenberg.de/download/TOGA/human_hg38_reference/"

    # Fetch all species
    species_directories = fetch_directory_listing(BASE_URL, args.animal_families)

    # Download alignment files
    all_downloaded_alignment_files = download_files(species_directories, 'proteinAlignments.fa.gz')

    # Filter long protein sequences
    filtered_files = filter_long_protein_sequences(all_downloaded_alignment_files, max_length=args.max_length)

    # Count the number of mapped orthologs for all species to humans
    all_mapped_ortholog_df = count_mapped_orthologs(filtered_files)

    # Remove proteins with more than max_orthologs in any species
    n_alignments = all_mapped_ortholog_df.value_counts(['protein', 'species']).unstack()
    n_alignments = n_alignments[(n_alignments > 0).all(axis=1)]
    proteins_to_keep = set(n_alignments.index[(n_alignments < args.max_orthologs).all(axis=1)])

    refiltered_files = filter_multimap_protein_sequences(filtered_files, proteins_to_keep)

    # Count the number of mapped orthologs after filtering
    filtered_mapped_ortholog_df = count_mapped_orthologs(refiltered_files)

    # Save the filtered DataFrame to a CSV file
    os.makedirs(args.output_folder, exist_ok=True)
    output_file = os.path.join(args.output_folder, "filtered_mapped_orthologs.csv")
    filtered_mapped_ortholog_df.to_csv(output_file, index=False)

    # Prepare dataset card data
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
        "data_file_path": "filtered_mapped_orthologs.csv"
    }

    # Create dataset card
    dataset_card = DatasetCard(args.template_path)

    # Check if user wants to publish the dataset
    if args.publish:
        # Get Hugging Face token from args or environment
        token = args.token or os.environ.get("HF_TOKEN")
        if not token:
            raise ValueError("Hugging Face token is required for publishing. Set it as an environment variable HF_TOKEN or pass it with --token.")

        # Push to Hugging Face Hub
        pusher = DatasetPusher(folder=args.output_folder, name="filtered-zoonomia-orthologs", token=token)
        pusher.push_to_hub(dataset_card, card_data)

        print("Filtered dataset creation and push to Hugging Face Hub completed.")
    else:
        print("Filtered dataset creation completed locally. No push to Hugging Face Hub.")

    # Clean up temporary files
    for file in all_downloaded_alignment_files + filtered_files + refiltered_files:
        os.remove(file)
    print("Temporary files cleaned up.")

if __name__ == "__main__":
    main()
