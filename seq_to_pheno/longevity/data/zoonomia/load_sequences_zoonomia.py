from fetch_sequences_zoonomia import SpeciesSequenceObject

if __name__ == "__main__":
    species = SpeciesSequenceObject.load('data', 'protein_sequence_df.tsv')
# load sequences