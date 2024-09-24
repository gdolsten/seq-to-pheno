---
language:
- en
license: {license}
size_categories:
- {size_category}
task_categories:
- {task_category}
pretty_name: {pretty_name}
dataset_info:
  features:
  - name: transcript
    dtype: string
  - name: protein
    dtype: string
  - name: mapped_to
    dtype: string
  - name: species
    dtype: string
  splits:
  - name: train
    num_bytes: {num_bytes}
    num_examples: {num_examples}
  download_size: {download_size}
  dataset_size: {dataset_size}
configs:
- config_name: default
  data_files:
  - split: train
    path: {data_file_path}
tags:
- zoonomia
- orthologs
---

# Zoonomia Orthologs Dataset

This dataset contains mapped orthologs for various species from the Zoonomia Project. It includes information about protein alignments and ortholog mappings between human and other species.

## Dataset Description

- **Point of Contact:** [Your Name or Organization]
- **Dataset Summary:** This dataset provides information about orthologous proteins across various species from the Zoonomia Project, mapped to human proteins.
- **Supported Tasks:** Comparative genomics, evolutionary studies, protein function prediction
- **Languages:** Not applicable (genomic data)
- **Source Data:** [Zoonomia Project](http://genome.senckenberg.de/download/TOGA/human_hg38_reference/)

## Dataset Structure

The dataset consists of a single CSV file with the following columns:

- `transcript`: Human transcript ID
- `protein`: Human protein name
- `mapped_to`: Non-human (query) organism transcript ID
- `species`: Name of the query species

## Data Splits

This dataset contains only a train split, which includes all the data.

## Dataset Creation

### Curation Rationale

This dataset was created to facilitate comparative genomics studies and to provide easy access to ortholog mappings across various species from the Zoonomia Project.

### Source Data

The source data was obtained from the Zoonomia Project, specifically from protein alignment files available at http://genome.senckenberg.de/download/TOGA/human_hg38_reference/.

### Annotations

The dataset does not contain additional annotations beyond the ortholog mappings provided in the original data.

## Considerations for Using the Data

### Social Impact of Dataset

This dataset can contribute to advancements in comparative genomics and evolutionary studies, potentially leading to better understanding of genetic diseases and evolutionary processes.

### Discussion of Biases

The dataset may have biases inherent in the original data collection and alignment processes used by the Zoonomia Project. Users should be aware of potential biases in species representation and alignment quality.

### Other Known Limitations

- The dataset only includes species available in the Zoonomia Project.
- Alignment quality may vary across different species and proteins.

## Additional Information

### Dataset Curators

[Your Name or Organization]

### Licensing Information

{license}

### Citation Information

If you use this dataset, please cite both this dataset and the original Zoonomia Project:

[Include citation for this dataset and the Zoonomia Project]

### Contributions

Thanks to the Zoonomia Project for making the original data available.