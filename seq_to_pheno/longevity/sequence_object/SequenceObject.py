import pandas as pd
import os
import numpy as np

class SpeciesSequenceObject():
    def load(file_path):
        SpeciesSequenceObject(pd.read_csv(file_path, sep='\t', index_col=[0, 1, 2]))


    def __init__(self, results_df):
        self.results_df = results_df
        self.all_organisms = self.results_df.index.get_level_values(0)
        self.all_protein_names = self.results_df.index.get_level_values(1)
    def select(self, *queries):
        subdf = self.results_df
        for query in queries:
            if query in self.all_organisms:
                subdf = subdf.xs(level = 0, key = query, axis=0, drop_level = False)
            elif query in self.all_protein_names:
                subdf = subdf.xs(level = 1, key = query, axis=0, drop_level = False)
        return subdf
    
    def saveas(self, basedir, file_path, shard=20):
        """Save the DataFrame into multiple shard files."""
        file_path = file_path.replace(".tsv", "")
        # Create shards of the DataFrame based on the shard size
        num_rows = len(self.results_df)
        shard_size = num_rows // shard if shard > 0 else num_rows
        
        # Ensure the directory exists for saving shards
        base_dir = file_path
        if base_dir and not os.path.exists(base_dir):
            os.makedirs(base_dir)
        
        # Save each shard to a separate file
        for i in range(shard):
            start_row = i * shard_size
            if i == shard - 1:  # Ensure the last shard captures any remainder
                subdf = self.results_df.iloc[start_row:]
            else:
                subdf = self.results_df.iloc[start_row:start_row + shard_size]
            
            shard_file_path = f"{basedir}/{file_path}/{file_path}_part_{i+1}.tsv"
            subdf.to_csv(shard_file_path, sep='\t', index=True)
            print(f"Shard {i+1} saved to {shard_file_path}")

    @classmethod
    def load(cls, basedir, file_path, shard=20):
        file_path = file_path.replace(".tsv", "")
        """Load multiple shard files and combine them into a single DataFrame."""
        shards = []
        
        # Load each shard and append to the list
        for i in range(shard):
            shard_file_path = f"{basedir}/{file_path}/{file_path}_part_{i+1}.tsv"
            if os.path.exists(shard_file_path):
                subdf = pd.read_csv(shard_file_path, sep='\t', index_col=[0, 1, 2])
                shards.append(subdf)
            else:
                print(f"Shard {i+1} not found: {shard_file_path}")
        
        # Concatenate all shards back into a single DataFrame
        if shards:
            full_df = pd.concat(shards)
            return cls(full_df)
        else:
            raise FileNotFoundError(f"No shard files found for base file: {file_path}")
        
    def to_fasta(self, basedir, file_path, protein_subset=None,):
        """
        Converts a set of sequences to FASTA format.

        Parameters:
        sequence_set (set): A set containing tuples of species and sequences.

        Returns:
        str: A string containing the sequences in FASTA format.
        """
        fastq = ""
        if protein_subset:
            data = self.select(protein_subset)
        else:
            data = self.results_df
        for _, row in data.reset_index().iterrows():
            fastq += f">{row.organism_name}.{row.protein_name}.{row.quer_transcript}\n{row.sequence}\n"

        with open(basedir + '/' + file_path, 'w') as file:
            file.write(fastq)
        return None
