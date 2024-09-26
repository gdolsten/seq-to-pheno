from sequence_object.SequenceObject import SpeciesSequenceObject

if __name__ == "__main__":
    species = SpeciesSequenceObject.load('data/zoonomia/', 'protein_sequence_df.tsv')
    # To get just TP53
    species.to_fasta('data/zoonomia/', 'TP53_protein_sequences.fasta', 'TP53')
    species.to_fasta('data/zoonomia/', 'all_protein_sequences.fasta')
    # To get all
    # species.to_fasta('data/zoonomia/', 'TP53_protein_sequences.fasta', 'TP53')
