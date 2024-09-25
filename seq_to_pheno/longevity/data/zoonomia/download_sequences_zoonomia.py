import requests
from bs4 import BeautifulSoup
from urllib.parse import urljoin
import os
import  gzip
import numpy as np

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



ANIMAL_FAMILIES = ['Afrotheria', 'Carnivora', 'Chiroptera', 'Dermoptera', 'Eulipotyphla', 'Lagomorpha', 'Metatheria', 
                'Perissodactyla', 'Pholidota', 'Primates', 'Prototheria', 'Rodentia', 'Ruminantia', 'Scandentia', 
                'Suina', 'Tylopoda', 'Whippomorpha', 'Xenarthra']

BASE_URL = "http://genome.senckenberg.de/download/TOGA/human_hg38_reference/"

# Fetch all the species
species_directories = fetch_directory_listing(BASE_URL, ANIMAL_FAMILIES)

# Download alignment files
all_downloaded_alignment_files = download_files(species_directories, 'proteinAlignments.fa.gz')

# Filter protein sequences longer than 1000 AA, remove them
filtered_files = filter_long_protein_sequences(all_downloaded_alignment_files, max_length=1000)

# Count the number of mapped orthologs for all species to humans
all_mapped_ortholog_df = count_mapped_orthologs(filtered_files)

MAX_NUMBER_ORTHOLOGS = 20
n_alignments = all_mapped_ortholog_df.value_counts(['protein', 'species']).unstack()
n_alignments = n_alignments[(n_alignments > 0).all(axis=1)]
proteins_to_keep = set(n_alignments.index[(n_alignments < MAX_NUMBER_ORTHOLOGS).all(axis=1)])

# Remove proteins with more than 20 orthologs in any species
refiltered_files = filter_multimap_protein_sequences(filtered_files, proteins_to_keep)



def fetch_sequences_from_files(file_paths, protein_name):
    """
    Fetches the reference and query sequences for a specific protein from a file.

    Parameters:
    file_path (str): The path to the file containing the sequences.
    protein_name (str): The middle name of the protein to search for.

    Returns:
    tuple: A tuple containing the reference and query sequences if found, else None.
    """
    sequences = set()
    for file_path in file_paths:
        print("Processing", file_path)
        species = file_path.split('/')[-2]
        # Open the file and read its contents
        with gzip.open(file_path, 'rt') as file:
            lines = file.readlines()
        
        # Iterate through the lines to find reference and query sequences for the given protein
        for i in range(len(lines)):
            reference_check = f"| PROT | REFERENCE" in lines[i] and protein_name in lines[i]
            query_check = f"| PROT | QUERY" in lines[i] and protein_name in lines[i]
            if reference_check:
                reference_sequence = lines[i + 1].strip().replace("-", '').replace('X', '')
                sequences.add(('Human', reference_sequence))
            elif query_check:
                query_sequence = lines[i + 1].strip().replace("-", '').replace('X', '')
                sequences.add((species, query_sequence))
    return sequences


def fetch_all_sequences_from_files(file_paths):
    """
    Fetches the reference and query sequences for a specific protein from a file.

    Parameters:
    file_path (str): The path to the file containing the sequences.
    protein_name (str): The middle name of the protein to search for.

    Returns:
    tuple: A tuple containing the reference and query sequences if found, else None.
    """
    seqdict = {}
    for file_path in file_paths:
        print("Processing", file_path)
        species = file_path.split('/')[-2]
        # Open the file and read its contents
        with gzip.open(file_path, 'rt') as file:
            lines = file.readlines()
        
        # Iterate through the lines to find reference and query sequences for the given protein
        for i in range(len(lines)):
            reference_check = f"| PROT | REFERENCE" in lines[i]
            query_check = f"| PROT | QUERY" in lines[i]
            if reference_check:
                reference_sequence = lines[i + 1].strip().replace("-", '').replace('X', '')
                seqdict[('Homo_sapiens', reference_sequence)] = lines[i].strip()
            elif query_check:
                query_sequence = lines[i + 1].strip().replace("-", '').replace('X', '')
                seqdict[(species, query_sequence)] = lines[i].strip()

    print("Processing seqdict")
    results = []
    for (species, sequence), metadata in seqdict.items():
        metadata_vals = metadata.split(" |")[0].split(".")
        if len(metadata_vals) > 3:
            continue
        else:   
            ref_transcript, protein_name, quer_transcript = metadata_vals
            ref_transcript = ref_transcript.replace(">", "")
            if species == 'Homo_sapiens':
                quer_transcript = ref_transcript
            results.append([species, ref_transcript, protein_name, quer_transcript, sequence])

    results_df = pd.DataFrame(results, columns = ['Species', 'ref_transcript', 'protein_name', 'quer_transcript', 'sequence'])
    return results_df
    # return seqdict #results_df

results_df = fetch_all_sequences_from_files(refiltered_files)
results_df['organism_name'] = results_df['Species'].apply(lambda x: x.split("__")[0] if x != 'Homo_sapiens' else x)
results_df = results_df.set_index(['organism_name', 'protein_name', 'quer_transcript'])

results_df = results_df[results_df['sequence'].apply(len) > 0]

obj = SpeciesSequenceObject(results_df)
obj.select("TP53")
obj.saveas('data', 'protein_sequence_df.tsv')





protein_name = 'TP53'
sequences = fetch_sequences_from_files(refiltered_files, protein_name)
to_fasta(sequences, 'examples/all_tp53.fa')



protein_name = 'DNAL1'
sequences = fetch_sequences_from_files(refiltered_files, protein_name)
to_fasta(sequences, 'examples/all_dnal1.fa')




protein_name = 'DNAL1'
proteins = ['ZNRF2', 'VDAC1']
import subprocess
for protein in proteins:
    sequences = fetch_sequences_from_files(refiltered_files, protein)
    to_fasta(sequences, f'examples/all_{protein}.fa')
    result = subprocess.run(['bash', 'run_esm.sh', protein], stdout=subprocess.PIPE)



def to_fasta(sequence_set, file_path):
    """
    Converts a set of sequences to FASTA format.

    Parameters:
    sequence_set (set): A set containing tuples of species and sequences.

    Returns:
    str: A string containing the sequences in FASTA format.
    """
    
    fastq = ""
    for c, (species, sequence) in enumerate(sequence_set):
        fastq += f">{species}.{c}\n{sequence}\n"

    with open(file_path, 'w') as file:
        file.write(fastq)
    
    return None


print((n_alignments > 0).sum(axis=1).sort_values(ascending=False))


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



