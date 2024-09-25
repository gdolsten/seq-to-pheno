import requests
from bs4 import BeautifulSoup
from urllib.parse import urljoin
import os
import  gzip
import pandas as pd
from huggingface_hub import HfApi
# from seq_to_pheno.hug.utils import DatasetCard, DatasetPusher
import sys
sys.path.append('C:/Users/MeMyself/seq-to-pheno/seq_to_pheno/hug')
from utils import DatasetCard, DatasetPusher
from typing import Optional

def fetch_directory_listing(url, animal_families):
    """Fetch and parse directory listing from a URL."""
    all_links = []
    for family in animal_families:
        family_url = url + family + "/"
        print(f"Fetching directory listing from {family_url}")
        response = requests.get(family_url, verify=False)
        soup = BeautifulSoup(response.text, 'html.parser')
        links = [urljoin(family_url, node.get('href')) for node in soup.find_all('a') if node.get('href').endswith('/')]
        links = [x for x in links if x != url] # Ignore base directory
        all_links += links
    return all_links

def download_files(urls, file_name='proteinAlignments.fa.gz'):
    """Download all files with a specific extension from the given directory URL."""
    all_downloaded_files = []
    for url in urls:
        file_url = url + file_name
        directory = url.split('/')[-2]  # Extract the directory name from the URL
        file_path = os.path.join('downloads', directory, file_name)  # Constructs the path where the file will be saved

        # Create directory if it doesn't exist
        os.makedirs(os.path.dirname(file_path), exist_ok=True,)

        # Start the download process
        if os.path.exists(file_path):
            print(f"File already exists: {file_path}")
        else:
            response = requests.get(file_url, stream=True, verify=False)
            if response.status_code == 200:
                print(file_path)
                with open(file_path, 'wb') as file:
                    for chunk in response.iter_content(chunk_size=8192):
                        file.write(chunk)
                print(f"Downloaded: {file_path}")
            else:
                print(f"Failed to download {file_url}: Status code {response.status_code}")
        all_downloaded_files.append(file_path)
    return all_downloaded_files

#  DEPRECATE: just work with gzipped files
# def decompress_gz_files(all_file_paths):
#     """Decompress all .gz files in the specified directory."""
#     decompressed_files = []
#     for filename in all_file_paths:
#             gz_path = filename
#             out_path = filename[:-3]
            
#             # Open the gzipped file and create a decompressed output file
#             with gzip.open(gz_path, 'rb') as f_in:
#                 with open(out_path, 'wb') as f_out:
#                     f_out.write(f_in.read())
#             print(f"Decompressed: {gz_path} -> {out_path}")
#             decompressed_files.append(out_path)
#     return decompressed_files


def count_mapped_orthologs(file_paths):
    """
    Counts the number of orthologs mapped to human for each species

    Parameters:
    file_path (str): The path to the file containing the sequences.

    Returns:
    pd.DataFrame: A DataFrame containing the number of orthologs mapped to human for each species.
        of the shape: [transcript, protein, mapped_to, species],
        where:
            transcript is the HUMAN (reference) transcript ID,
            protein is the HUMAN (reference) protein name,
            mapped_to is non-human (query) organism transcript ID (e.g. elephant's transcript ID)
            species is the name of the query species (e.g. 'Elephas_maximus__Asiatic_elephant__HLeleMax1'),
            
            The values in transcript and mapped_to are the orthologs.

        For example, the first row should be:
            transcript          protein        mapped_to   species
            >ENST00000409217    ZFAND2B        77          Chrysochloris_asiatica__Cape_golden_mole__chrAsi1
            
            In this case, ENST00000409217<->77, so ENST00000409217 is the 
                human transcript ID which is mapped to transcript 77 in the Cape golden mole.
            ENST00000409217 codes for ZFAND2B
        Note: 
            - The file contains the mapping of orthologs to human, so the human transcript ID is always the reference.
            - The file contains the mapping of proteins, not genes.
            - A protein can be mapped to multiple orthologs in the same species.
            - A protein can be mapped to multiple transcripts in humans
    """
    all_mapped_orthologs = []
    for file_path in file_paths:
        species = file_path.split('/')[-2]
        mapped_orthologs = []
        print("Processing", species)
        with gzip.open(file_path, 'rt') as infile:
            lines = infile.readlines()

            # Iterate through the lines to find reference and query sequences
            i = 0
            while i < len(lines):
                assert lines[i].startswith('>') and ("| PROT | REFERENCE" in lines[i] or "| PROT | QUERY" in lines[i])
                header1 = lines[i].strip()
                metadata = header1.split('|')[0].strip().split(".")
                if len(metadata) > 3:
                    print("Skipping", metadata)
                    pass
                else:
                    transcript, protein, mapped_to = metadata
                    mapped_orthologs.append([transcript, protein, mapped_to,])
                i += 4  # Skip to the next sequence block
        mapped_orthologs = pd.DataFrame(mapped_orthologs, columns=['transcript', 'protein', 'mapped_to'])
        mapped_orthologs['species'] = species
        all_mapped_orthologs.append(mapped_orthologs)
    all_mapped_orthologs = pd.concat(all_mapped_orthologs, axis=0)
    return all_mapped_orthologs

def filter_long_protein_sequences(file_paths, max_length=1000):
    """
    Filters sequences in the input file to only include those with fewer than `max_length` amino acids,
    and writes the filtered sequences to a new file.

    Parameters:
    file_path (str): The path to the file containing the sequences.
    max_length (int): The maximum allowed length for sequences (default is 1000).
    """
    output_files = []
    for file_path in file_paths:
        print("Processing", file_path)
        output_file = os.path.dirname(file_path) + f"/filtered.maxlen={max_length}." + os.path.basename(file_path)
        if os.path.exists(file_path):
            output_files.append(output_file)
            # if os.path.exists(output_file):
                # print(f"File already exists: {output_file}")
                # continue
            with gzip.open(file_path, 'rt') as infile, gzip.open(output_file, 'wt') as outfile:
                lines = infile.readlines()

                # Iterate through the lines to find reference and query sequences
                i = 0
                while i < len(lines)-6:
                    if lines[i].startswith('>') and ("| PROT | REFERENCE" in lines[i] or "| PROT | QUERY" in lines[i]):
                        sequence1 = lines[i + 1].strip()
                        sequence2 = lines[i + 3].strip()
                        # Check if the sequence length is less than or equal to max_length
                        if (len(sequence1) <= max_length) or (len(sequence2) <= max_length):
                            # Write the header and sequence to the output file
                            outfile.write(lines[i])  # Write the header
                            outfile.write(lines[i+1])  # Write the sequence
                            outfile.write(lines[i+2])  # Write the header
                            outfile.write(lines[i+3])  # Write the sequence
                        else:
                            pass
                    i += 4  # Skip to the next sequence block
    return output_files

def filter_multimap_protein_sequences(file_paths, proteins_to_keep):
    """
    Filters sequences in the input file to only include those with fewer than `max_number` orthologs,
    and writes the filtered sequences to a new file.

    Parameters:
    file_path (str): The path to the file containing the sequences.
    max_length (int): The maximum allowed length for sequences (default is 1000).
    """
    output_files = []
    for file_path in file_paths:
        print("Processing", file_path)
        success = 0
        failure = 0
        output_file = os.path.dirname(file_path) + f"/nomultimap." + os.path.basename(file_path)
        output_files.append(output_file)
        with gzip.open(file_path, 'rt') as infile, gzip.open(output_file, 'wt') as outfile:
            lines = infile.readlines()

            # Iterate through the lines to find reference and query sequences
            i = 0
            while i < len(lines):
                if lines[i].startswith('>') and ("| PROT | REFERENCE" in lines[i] or "| PROT | QUERY" in lines[i]):
                    # Check if the protein has few enough orthologs
                    protein = lines[i].strip().split("|")[0].split('.')[1]
                    if (protein in proteins_to_keep):
                        # Write the header and sequence to the output file
                        outfile.write(lines[i])  # Write the header
                        outfile.write(lines[i+1])  # Write the sequence
                        outfile.write(lines[i+2])  # Write the header
                        outfile.write(lines[i+3])  # Write the sequence
                        success += 1
                    else:
                        failure += 1
                    
                i += 4  # Skip to the next sequence block
        print(f"Success: {success}, Failure: {failure}; {success/(success+failure):.2f} success rate")
    return output_files

# refiltered_files = filter_multimap_protein_sequences(filtered_files, proteins_to_keep)


# # ANIMAL_FAMILIES = ['Afrotheria', 'Carnivora', 'Chiroptera', 'Dermoptera', 'Eulipotyphla', 'Lagomorpha', 'Metatheria', 
# #                 'Perissodactyla', 'Pholidota', 'Primates', 'Prototheria', 'Rodentia', 'Ruminantia', 'Scandentia', 
# #                 'Suina', 'Tylopoda', 'Whippomorpha', 'Xenarthra']

# # BASE_URL = "http://genome.senckenberg.de/download/TOGA/human_hg38_reference/"

# # # Fetch all the species
# # species_directories = fetch_directory_listing(BASE_URL, ANIMAL_FAMILIES)

# # # Download alignment files
# # all_downloaded_alignment_files = download_files(species_directories, 'proteinAlignments.fa.gz')

# # Filter protein sequences longer than 1000 AA, remove them
# filtered_files = filter_long_protein_sequences(all_downloaded_alignment_files, max_length=1000)

# # Count the number of mapped orthologs for all species to humans
# all_mapped_ortholog_df = count_mapped_orthologs(filtered_files)

# MAX_NUMBER_ORTHOLOGS = 20
# n_alignments = all_mapped_ortholog_df.value_counts(['protein', 'species']).unstack()
# n_alignments = n_alignments[(n_alignments > 0).all(axis=1)]
# proteins_to_keep = set(n_alignments.index[(n_alignments < MAX_NUMBER_ORTHOLOGS).all(axis=1)])
# # Remove proteins with more than 20 orthologs in any species
# refiltered_files = filter_multimap_protein_sequences(filtered_files, proteins_to_keep)



# # Filter long protein sequences, remove them
# filtered_files = filter_long_protein_sequences(all_downloaded_alignment_files, max_length=1000)

# # Count the number of mapped orthologs for all species to humans
# all_mapped_ortholog_df = count_mapped_orthologs(filtered_files, max_length=1000)

# n_alignments = all_mapped_ortholog_df.value_counts(['protein', 'species']).unstack()

# print((n_alignments > 0).sum(axis=1).sort_values(ascending=False))


# import pandas as pd
# import numpy as np
# data = pd.read_csv('downloads/Elephas_maximus__Asiatic_elephant__HLeleMax1/proteinAlignments.fa', sep='\t')


# # Fetch BRCA protein sequences and print them
# # brca_sequences, fasta_data = fetch_brca_sequences()


# filter_long_protein_sequences('downloads/Elephas_maximus__Asiatic_elephant__HLeleMax1/proteinAlignments.fa', max_length=1000)



# import gzip
# import os
                
# def fetch_lengths_from_file(file_path, protein_name):
#     """
#     Fetches the reference and query sequences for a specific protein from a file.

#     Parameters:
#     file_path (str): The path to the file containing the sequences.
#     protein_name (str): The middle name of the protein to search for.

#     Returns:
#     tuple: A tuple containing the reference and query sequences if found, else None.
#     """
#     reference_sequence = ""
#     query_sequence = ""
    
#     # Open the file and read its contents
#     with open(file_path, 'r') as file:
#         lines = file.readlines()
    
#     # Iterate through the lines to find reference and query sequences for the given protein
#     seqlens = []
#     for i in range(len(lines)):
#         if '>' in lines[i]:
#             pass
#         else:
#             seqlens.append(len(lines[i]))
#     return seqlens



def fetch_sequences_from_file(file_path, protein_name):
    """
    Fetches the reference and query sequences for a specific protein from a file.

    Parameters:
    file_path (str): The path to the file containing the sequences.
    protein_name (str): The middle name of the protein to search for.

    Returns:
    tuple: A tuple containing the reference and query sequences if found, else None.
    """
    reference_sequence = ""
    query_sequence = ""
    
    # Open the file and read its contents
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    # Iterate through the lines to find reference and query sequences for the given protein
    for i in range(len(lines)):
        reference_check = f"| PROT | REFERENCE" in lines[i] and protein_name in lines[i]
        if reference_check:
            if len(reference_sequence) > 0:
                raise ValueError("Multiple reference sequences found for the same protein.")
            reference_sequence = lines[i + 1].strip()
        elif f"| PROT | QUERY" in lines[i] and protein_name in lines[i]:
            if len(query_sequence) > 0:
                raise ValueError("Multiple reference sequences found for the same protein.")            
            query_sequence = lines[i + 1].strip()
    
    if reference_sequence and query_sequence:
        return reference_sequence, query_sequence
    else:
        return None


# Hugging Face dataset creation class
class DatasetPusher:
    def __init__(self, folder: str, name: str, token: Optional[str] = None):
        self.folder = folder
        self.name = name
        self.token = token or os.environ.get("HF_TOKEN")
        if not self.token:
            raise ValueError("Hugging Face token is required. Set it as an environment variable HF_TOKEN or pass it to the constructor.")
        self.api = HfApi()

    def push_to_hub(self):
        current_time = datetime.now().strftime("%Y%m%d_%H%M%S")
        repo_name = f"{self.name}-{current_time}"

        try:
            self.api.repo_info(repo_id=repo_name, repo_type="dataset")
            print(f"Repository {repo_name} already exists. Pushing to existing repository.")
        except Exception:
            self.api.create_repo(repo_id=repo_name, repo_type="dataset", private=False)
            print(f"Created new repository: {repo_name}")

        local_dir = f"./temp_{repo_name}"
        repo = Repository(local_dir=local_dir, clone_from=f"{self.api.whoami()['name']}/{repo_name}", repo_type="dataset", use_auth_token=self.token)

        for root, _, files in os.walk(self.folder):
            for file in files:
                src_path = os.path.join(root, file)
                dst_path = os.path.join(local_dir, os.path.relpath(src_path, self.folder))
                os.makedirs(os.path.dirname(dst_path), exist_ok=True)
                with open(src_path, "rb") as src, open(dst_path, "wb") as dst:
                    dst.write(src.read())

        repo.git_add()
        repo.git_commit("Initial commit")
        repo.git_push()

        print(f"Successfully pushed dataset to https://huggingface.co/datasets/{self.api.whoami()['name']}/{repo_name}")

        repo.delete_repo()
        print("Cleaned up local repository.")

if __name__ == "__main__":
    ANIMAL_FAMILIES = ['Afrotheria', 'Carnivora', 'Chiroptera', 'Dermoptera', 'Eulipotyphla', 'Lagomorpha', 'Metatheria', 
                    'Perissodactyla', 'Pholidota', 'Primates', 'Prototheria', 'Rodentia', 'Ruminantia', 'Scandentia', 
                    'Suina', 'Tylopoda', 'Whippomorpha', 'Xenarthra']

    BASE_URL = "http://genome.senckenberg.de/download/TOGA/human_hg38_reference/"

    # Fetch all the species
    species_directories = fetch_directory_listing(BASE_URL, ANIMAL_FAMILIES)

    # Download alignment files
    all_downloaded_alignment_files = download_files(species_directories, 'proteinAlignments.fa.gz')

    # Filter long protein sequences, remove them
    filtered_files = filter_long_protein_sequences(all_downloaded_alignment_files, max_length=1000)

    # Count the number of mapped orthologs for all species to humans
    all_mapped_ortholog_df = count_mapped_orthologs(filtered_files)

    # Remove proteins with more than 20 orthologs in any species
    MAX_NUMBER_ORTHOLOGS = 20

    
    n_alignments = all_mapped_ortholog_df.value_counts(['protein', 'species']).unstack()
    n_alignments = n_alignments[(n_alignments > 0).all(axis=1)]
    proteins_to_keep = set(n_alignments.index[(n_alignments < MAX_NUMBER_ORTHOLOGS).all(axis=1)])

    refiltered_files = filter_multimap_protein_sequences(filtered_files, proteins_to_keep)

    # Count the number of mapped orthologs for all species to humans after filtering
    filtered_mapped_ortholog_df = count_mapped_orthologs(refiltered_files)
    # Save the filtered DataFrame to a CSV file
    output_folder = "filtered_zoonomia_dataset"
    os.makedirs(output_folder, exist_ok=True)
    output_file = os.path.join(output_folder, "filtered_mapped_orthologs.csv")
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
        "data_file_path": "mapped_orthologs.csv"
    }

    # Create dataset card
    template_path = "seq_to_pheno/hug/zoonomia_dataset_repo_template/README.md"
    dataset_card = DatasetCard(template_path)
    # Update the dataset card with filtering information
    # dataset_card.add_section("Filtering Details", f"""
    # This dataset has been filtered to remove:
    # 1. Proteins longer than 1000 amino acids
    # 2. Proteins with more than {MAX_NUMBER_ORTHOLOGS} orthologs in any species

    # Original number of mapped orthologs: {len(all_mapped_ortholog_df)}
    # Filtered number of mapped orthologs: {num_examples}
    # """)
    #token=os.environ.get("HF_TOKEN")
    token = "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"
    # Create and push the dataset to Hugging Face Hub
    # Create and push the dataset to Hugging Face Hub
    # token = os.environ.get("HF_TOKEN")
    if not token:
        raise ValueError("Hugging Face token is required. Set it as an environment variable HF_TOKEN.")

    pusher = DatasetPusher(folder=output_folder, name="filtered-zoonomia-orthologs", token=token)
    pusher.push_to_hub(dataset_card, card_data)

    print("Filtered dataset creation and push to Hugging Face Hub completed.")

    # Clean up temporary files
    for file in all_downloaded_alignment_files + filtered_files + refiltered_files:
        os.remove(file)
    print("Temporary files cleaned up.")
    