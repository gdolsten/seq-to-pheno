#!/bin/bash

# Assign the first command line argument to the variable 'protein_name'
protein_name=$1

# Use the variable in the file paths
esm-extract esm2_t33_650M_UR50D examples/all_${protein_name}.fa examples/test_${protein_name}_emb --repr_layers 0 32 33 --include mean
