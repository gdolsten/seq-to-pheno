# seq-to-pheno

Code for longevity : Put in longevity project

Code for depmap : Put in depmap project

Fetching embeddings from sequence : get_embeddings.py

## Datasets

- [add datasets here](https://huggingface.co/seq-to-pheno)

### Get

**get the filtered ortholog dataset :**

```sh
curl -X GET \
     "https://datasets-server.huggingface.co/first-rows?dataset=seq-to-pheno%2Ffiltered_orthologs&config=default&split=train"
```

**get the mapped ortholog dataset :**

```sh 
curl -X GET \
     -H "Authorization: Bearer $HF_TOKEN" \
     "https://datasets-server.huggingface.co/rows?dataset=seq-to-pheno%2Fmapped_orthologs&config=default&split=train&offset=0&length=100"
```

### Use

```python
from datasets import load_dataset

mapped = load_dataset("seq-to-pheno/mapped_orthologs")
```

```python
from datasets import load_dataset

mapped = load_dataset("seq-to-pheno/filtered_orthologs")
```


### Re-Create the filtered Ortholog Dataset:

```sh
python ./scripts/filtered_dataset.py --folder /downloads --template_path /seq_to_pheno/hug/zoonomia_dataset_repo_template/README.md --token hf_xxx --max_length 1000 --max_orthologs 20 --publish
```

### Re-Create the Fasta Zoonotica Dataset:

**To extract sequences for a specific gene and publish:**

```sh
python extract_and_publish_protein_sequences.py --input_folder data/zoonomia/ --input_file protein_sequence_df.tsv --output_folder data/zoonomia/ --output_file TP53_protein_sequences.fasta --gene TP53 --publish --repo_name filtered-zoonomia-tp53 --hf_token hf_your_token
```

**To extract all sequences and publish:**

```sh
python extract_and_publish_protein_sequences.py --input_folder data/zoonomia/ --input_file protein_sequence_df.tsv --output_folder data/zoonomia/ --output_file all_protein_sequences.fasta --publish --repo_name filtered-zoonomia-all --hf_token hf_your_token
```
