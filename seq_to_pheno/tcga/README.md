# ICGC/TCGA Preprocessing

Space for processing ICGC/TCGA data

### Source
The ICGC (International Cancer Genome Consortium) in collaboration with TCGA and the Sanger Institute has created a hub for a wide-range of NGS data and metadata across various cancer types.

Instructions on how to access the data are [found on this page](https://docs.icgc-argo.org/docs/data-access/icgc-25k-data).

## Requirements
#### For you Shell environment
```
wget
unzip
java (openjdk-23)
tabix
bcftools
samtools
git
```

#### Python packages:
```
collections
pandas
pysam
intervaltree
biopython
requests
tqdm
```

## ETL Overview
**get normalized transcript expression data**

```
aws s3 ls s3://icgc25k-open/PCAWG/transcriptome/transcript_expression/pcawg.rnaseq.transcript.expr.tpm.tsv.gz --endpoint-url https://object.genomeinformatics.org --no-sign-request
```

This is transcriptomics data that has been normalized in a manner (trascripts per million) that make it comparable between subjects -- which is great for use in models comparing outcomes between subjects of with varied genotypes/phenotypes.

**Table Structure**:
- Column 1: Transcript ID
- Column 2+: Metadata-linked aliquot_id


**get aliquot_id linked metadata**
```
aws s3 ls s3://icgc25k-open/PCAWG/transcriptome/metadata/rnaseq.extended.metadata.aliquot_id.V4.tsv.gz --endpoint-url https://object.genomeinformatics.org --no-sign-request
```

This table links the aliquot_id to a plethora of information about the sample.

**Sample headers:**
```
aliquot_id      study   project_code    icgc_donor_id   icgc_specimen_id        icgc_sample_id  submitted_donor_id      submitted_specimen_id   submitted_sample_id     specimen_type   donor_id        specimen_id     sample_id       date    file_name filesize        md5_checksum    analysis_id     uri     center_name     platform        platform_model  lib_id  rg_label        files   files_blacklist pilot_63        includes_spike_ins      spike_ins_fasta spike_ins_concentration library_type      library_strategy        accession       state   is_tumour       file_type       unique_donor_id wgs_white_black_gray    wgs_exclusion_white_gray        donor_unique_id tcga_donor_uuid donor_sex       donor_vital_status      donor_diagnosis_icd10     first_therapy_type      first_therapy_response  donor_age_at_diagnosis  donor_survival_time     donor_interval_of_last_followup tobacco_smoking_history_indicator       tobacco_smoking_intensity       alcohol_history alcohol_history_intensity donor_wgs_included_excluded     tcga_specimen_uuid      tcga_sample_uuid        organ_system    histology_abbreviation  histology_tier1 histology_tier2 histology_tier3 histology_tier4 tumour_histological_code        tumour_histological_type  tumour_stage    tumour_grade    percentage_cellularity  level_of_cellularity    tcga_expert_re.review   tumour_histological_comment     specimen_donor_treatment_type   rna_seq_specimen_id     Color (RGB code)        rna_seq_aliquot_id        wgs_aliquot_id  GTex.tissue     tumor.normal    matched_wgs_aliquot_id
```

When thinking about linking transcriptomics data to longevity we have many fields here to choose from:
- aliquot_id (matched Id from metadata to transcript data)
- donor_survival_time (days)
- project_code (ICGC Cancer Symbol)
- GTex.tissue (Layman Cancer Type)

**Example:**
```
➜ head -5 rnaseq.extended.metadata.aliquot_id.V4.tsv | cut -f1,3,75,76,48 
aliquot_id      project_code    donor_survival_time     wgs_aliquot_id  GTex.tissue
b337121c-9821-4644-820e-b8c477f6c38a    GBM-US  489     d60f54f5-b154-42c4-99fb-cea4e7a33dc7    Brain
612ef912-5a28-4c11-8703-3376f51afef5    GBM-US  317     c065761d-f775-457f-bda0-4c7c257a701e    Brain
56a705f4-fd28-44ff-8a3c-85bc4300c760    GBM-US  737     56ffaa35-814c-4c0b-b3c6-d4514d34fec2    Brain
28239a0e-9990-49ef-a159-ba63fb078c77    GBM-US  485     7cae6c0b-36fe-411b-bbba-093a4c846d84    Brain
```

**Get Linked SNV and Indel Data**

The ICGC data collection also provides VCF files for indels and SNVs (single nucleotide variants)

```
aws s3 cp s3://icgc25k-open/PCAWG/consensus_snv_indel/final_consensus_snv_indel_passonly_icgc.public.tgz --endpoint-url https://object.genomeinformatics.org --no-sign-request ./

tar -zxvf final_consensus_snv_indel_passonly_icgc.public.tgz
```

The extracted tar.gz file gives us two separate directories for SNVs and Indels:

```
indel/
snv_mnv/
```

Here is a sneak-peak inside the snv_mv directory:

```
➜ ls snv_mnv | head
0009b464-b376-4fbc-8a56-da538269a02f.consensus.20160830.somatic.snv_mnv.vcf.gz 
0009b464-b376-4fbc-8a56-da538269a02f.consensus.20160830.somatic.snv_mnv.vcf.gz.tbi 
003819bc-c415-4e76-887c-931d60ed39e7.consensus.20160830.somatic.snv_mnv.vcf.gz 
003819bc-c415-4e76-887c-931d60ed39e7.consensus.20160830.somatic.snv_mnv.vcf.gz.tbi 
```

The prefix for each of these files corresponds to the wgs_aliquot_id from the metadata file. Meaning we can link our transcripts data to genotype data.

And here's a look at the format of one of these VCF files (after decompression):

```
#CHROM  POS ID  REF ALT QUAL    FILTER  INFO
1   1002242 .   C   T   .   .   Callers=broad,dkfz,muse,sanger;NumCallers=4;VAF=0.5;cosmic=COSN8882299;t_alt_count=14;t_ref_count=14;Variant_Classification=RNA
1   4400126 .   C   T   .   .   1000genomes_AF=0.00019968;1000genomes_ID=rs188587967;Callers=broad,dkfz,muse,sanger;NumCallers=4;VAF=0.2;cosmic=COSN8470412;dbsnp=rs188587967;t_alt_count=6;t_ref_count=24;Variant_Classification=IGR
1   6544147 .   C   T   .   .   Callers=broad,dkfz,muse,sanger;NumCallers=4;VAF=0.3514;cosmic=COSN8504060;t_alt_count=13;t_ref_count=24;Variant_Classification=Intron
1   6999251 .   G   A   .   .   Callers=broad,dkfz,muse,sanger;NumCallers=4;VAF=0.5455;cosmic=COSN8505688;t_alt_count=18;t_ref_count=15;Variant_Classification=Intron
1   7043470 .   C   G   .   .   Callers=broad,dkfz,muse,sanger;NumCallers=4;VAF=0.4643;cosmic=COSN8505922;repeat_masker=L2c;t_alt_count=13;t_ref_count=15;Variant_Classification=Intron
1   7947718 .   C   T   .   .   Callers=broad,dkfz,muse,sanger;NumCallers=4;VAF=0.9048;cosmic=COSN8510046;repeat_masker=AluJb;t_alt_count=19;t_ref_count=2;Variant_Classification=IGR
1   9585830 .   G   A   .   .   Callers=broad,dkfz,muse,sanger;NumCallers=4;VAF=0.5172;cosmic=COSN8522351;repeat_masker=LTR2;t_alt_count=15;t_ref_count=14;Variant_Classification=IGR
```

What we see here is information about the location of variants found in this subject, along with relevant information like alt_count and ref_count (which can be used to define allele depth and genotypes), as well as information related to the Variant_Classification, and the type of sequencing used to generate the data.

### Annotate Variants
For downstream analysis we will need to filter our variants based on their predicted protein-level effect. To do so, we will annotate our variants using SnpEff for variant predictions.

This will require some extra steps to prepare. We must make sure that we are using the correct reference genome (as described in the PCAWG documentation). We also need to make sure we are using the correct reference GTF -- which will be needed to construct the SnpEff database. 

And then finally, we will need to make sure to add contig-tags to each VCF file. The contig-tags are created from the reference genome used to make the original variant calls (GRCh37.75).

All of the steps above can be run using ```set_tcga_data.py``` **just make sure to adjust your own filepaths**.

>NOTE: You will need access to both bcftools and samtools. If you have any issues setting up -- feel free to reach out to me (Harrison/Almuraqib) or tinker on your own.

> NOTE: This script will take a long time to run, it might make sense to subset the list to only a handful of subjects first for downstream testing purposes.

### Merge Annotations
Once we have successfully preprocessed our annotations, we will need to combine the indels and snvs into a single vcf file. 

We can run this using ```merge_annotations.py```

### Get sequences
We can then take these filtered and merged variant calls and then generate protein sequences based on the variants found. It will account for all overlapping variants in a transcript and create a protein-fasta file of the resulting protein-sequence. 

We can run this using ```get_sequences.py```