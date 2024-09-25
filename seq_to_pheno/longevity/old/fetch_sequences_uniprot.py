import requests

def fetch_brca_sequences():
    # Define the API endpoint with the query for BRCA proteins in any species
    url = 'https://rest.uniprot.org/uniprotkb/search?query=brca1&format=fasta'
    
    try:
        # Send a GET request to the API
        response = requests.get(url)
        # Raise an error for a failed request
        response.raise_for_status()
        
        # Get the FASTA sequences as text
        fasta_data = response.text
        return parse_fasta(fasta_data), fasta_data
    except requests.exceptions.RequestException as e:
        print(f"An error occurred: {e}")
        return None


## has some hacks
def parse_fasta_header(header):
    """
    Parses a FASTA header and returns the accession number, description, and species.
    
    Args:
        header (str): The FASTA header line.
        
    Returns:
        dict: A dictionary with 'accession', 'description', and 'species' keys.
    """
    # Example header: >sp|O75417|DPOLQ_HUMAN DNA polymerase theta OS=Homo sapiens OX=9606 GN=POLQ PE=1 SV=2
    if not header.startswith(">"):
        raise ValueError("Invalid FASTA header format")

    # Remove the initial '>' character
    header = header[1:]

    # Split the header into parts using pipe (|) to separate database, accession, and entry name
    parts = header.split("|")
    
    if len(parts) < 3:
        raise ValueError("Unexpected FASTA header format")
    
    # Get the accession and description
    accession = parts[1]  # e.g., O75417
    description = parts[2].split(" OS=")[0]  # e.g., DPOLQ_HUMAN DNA polymerase theta
    
    # Extract the species (OS stands for Organism Species)
    species_info = header.split("OS=")[1].split("OX=")[0].strip()  # e.g., Homo sapiens
    
    # Extract the gene name (GN stands for Gene Name)
    gene_name = None
    if "GN=" in header:
        gene_name = header.split("GN=")[1].split()[0]  # e.g., POLQ
    
    return {
        "accession": accession,
        "description": description,
        "species": species_info,
        "gene_name": gene_name
    }

def parse_fasta(fasta_data):
    """
    Parses FASTA formatted data and returns a dictionary of accession numbers and sequences.

    Args:
        fasta_data (str): The raw FASTA data.

    Returns:
        dict: A dictionary with protein accession numbers as keys and sequences as values.
    """
    sequences = {}
    gene_key = None
    current_sequence = []
    
    for line in fasta_data.splitlines():
        if line.startswith(">"):
            if gene_key and current_sequence:
                sequences[gene_key] = ''.join(current_sequence)
            header_info = parse_fasta_header(line)
            accession = header_info["accession"]
            species = header_info["species"]
            gene_name = header_info["gene_name"]
            gene_key = (accession, species, gene_name)
            current_sequence = []
        else:
            current_sequence.append(line.strip())
    
    if gene_key and current_sequence:
        assert gene_key not in sequences
        sequences[gene_key] = ''.join(current_sequence)
    else:
        print("No sequences found in the provided data")

    return sequences


# Fetch BRCA protein sequences and print them
brca_sequences, fasta_data = fetch_brca_sequences()


if brca_sequences:
    print(brca_sequences)


for i in brca_sequences.split('\n')[:3]:
    print(i)