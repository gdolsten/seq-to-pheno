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

Once finished we can now see the **ANN** tag now added to the info field on our VCF files:

In the header:
```
##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'">
```

Example VCF line:
```
#CHROM  POS ID  REF ALT QUAL    FILTER  INFO
1   1273949 .   A   G   .   .   Callers=broad,dkfz,muse,sanger;NumCallers=4;VAF=0.3182;cosmic=COSM4909256;t_alt_count=28;t_ref_count=60;Variant_Classification=Missense_Mutation;ANN=G|missense_variant|MODERATE|DVL1|ENSG00000107404|transcript|ENST00000378888|protein_coding|12/15|c.1292T>C|p.Ile431Thr|1577/3239|1292/2088|431/695||,G|missense_variant|MODERATE|DVL1|ENSG00000107404|transcript|ENST00000378891|protein_coding|12/15|c.1217T>C|p.Ile406Thr|1264/2926|1217/2013|406/670||,G|downstream_gene_variant|MODIFIER|TAS1R3|ENSG00000169962|transcript|ENST00000339381|protein_coding||c.*4105A>G|||||3263|,G|downstream_gene_variant|MODIFIER|DVL1|ENSG00000107404|transcript|ENST00000472445|retained_intron||n.*3201T>C|||||3201|
```


### Merge Annotations
Once we have successfully preprocessed our annotations, we will need to combine the indels and snvs into a single vcf file. 

We can run this using ```merge_annotations.py```

This will create the directory ```seq_to_pheno/tcga/data/variants/vcf/``` that will contain our combined annotate vcf.

### Get sequences
We can now take these filtered and merged variant calls and then generate protein sequences based on the variants found. It will account for all overlapping variants in a transcript and create a protein-fasta file of the resulting protein-sequence. 

We can run this using ```get_sequences.py```

This will create two new directories:

```
seq_to_pheno/tcga/data/variants/mutated_proteins
seq_to_pheno/tcga/data/variants/wildtype_proteins/
```

These are dirrectories that will contain the various protein-fasta sequences of our mutated proteins as well as their wildtype counterparts.

Mutant transcripts
```
➜ ls seq_to_pheno/tcga/data/variants/mutated_proteins | head
0009b464-b376-4fbc-8a56-da538269a02f_ENST00000003084_mutated.fasta 
0009b464-b376-4fbc-8a56-da538269a02f_ENST00000020673_mutated.fasta 
0009b464-b376-4fbc-8a56-da538269a02f_ENST00000046087_mutated.fasta 
0009b464-b376-4fbc-8a56-da538269a02f_ENST00000192314_mutated.fasta 
0009b464-b376-4fbc-8a56-da538269a02f_ENST00000228327_mutated.fasta 
0009b464-b376-4fbc-8a56-da538269a02f_ENST00000229030_mutated.fasta 
0009b464-b376-4fbc-8a56-da538269a02f_ENST00000230859_mutated.fasta 
0009b464-b376-4fbc-8a56-da538269a02f_ENST00000245255_mutated.fasta 
0009b464-b376-4fbc-8a56-da538269a02f_ENST00000248058_mutated.fasta 
0009b464-b376-4fbc-8a56-da538269a02f_ENST00000253001_mutated.fasta

```

In mutated_proteins you will also find log files for variants that passed filtration but coulnd't be included in the protein transcript for one reason or another:

```
➜ ls seq_to_pheno/tcga/data/variants/mutated_proteins/*failed* | head
0009b464-b376-4fbc-8a56-da538269a02f_failed_variants.log 
005794f1-5a87-45b5-9811-83ddf6924568_failed_variants.log 
00b9d0e6-69dc-4345-bffd-ce32880c8eef_failed_variants.log 
00db1b95-8ca3-4cc4-bb46-6b8c8019a7c7_failed_variants.log 
0168a2a6-c3af-4d58-a51c-d33f0fc7876d_failed_variants.log 
02917220-6a7a-46a1-8656-907e96bef88e_failed_variants.log 
03ad38a6-0902-4aaa-84a3-91ea88fa9883_failed_variants.log 
03c3c692-8a86-4843-85ae-e045f0fa6f88_failed_variants.log 
046d7386-95c8-4501-9e55-c85bec272a7a_failed_variants.log 
04b570c2-3224-4e9b-81cc-089b4a7ff07a_failed_variants.log
```

>**NOTE**: Although these errors represent a small minority of all variant calls we process here -- these instances will need to be investigated eventually. Most likely, there is some mismatch with the reference used to create this call vs the one we've used so far. This is likely due to an upstream issue inherenet to the dataset, but that will have to be proven.


Wildtype transcripts:

```
➜ ls seq_to_pheno/tcga/data/variants/wildtype_proteins | head 
ENST00000001008.fasta 
ENST00000001146.fasta 
ENST00000002125.fasta 
ENST00000002165.fasta 
ENST00000002596.fasta 
ENST00000002829.fasta 
ENST00000003084.fasta 
ENST00000003302.fasta 
ENST00000003583.fasta 
ENST00000004103.fasta

```

### Protein Metadata Construction

Okay, now that we have all of our sequences generated we must them back to sample metadata. 

For this, I have created ```create_protein_metadata_table.py```

This script will create a table that will create a unique row for each mutated_protein we've generated and with columns that contain the related wildtype_protein as well as important metadata related to the indivdual that whoose tumor transcribed the mutant protein. 

This script can be run to only include relative paths to the proteins:

```
➜ python seq-to-pheno/seq_to_pheno/tcga/create_protein_metadata_table.py

Metadata table saved to seq_to_pheno/tcga/data/protein_sequences_metadata.tsv
```

Inside table:
```
aliquot_id      mutated_protein_path    wildtype_protein_path   wgs_aliquot_id  Cancer Type     Cancer Stage    Donor Survival Time     Donor Vital Status      Donor Age at Diagnosis Tumour Grade    Donor Sex       Histology Abbreviation
8fb9496e-ddb8-11e4-ad8f-5ed8e2d07381    seq_to_pheno/tcga/data/variants/mutated_proteins/80ab6c08-c622-11e3-bf01-24c6515278c0_ENST00000512632_mutated.fasta     seq_to_pheno/tcga/data/variants/wildtype_proteins/ENST00000512632.fasta        80ab6c08-c622-11e3-bf01-24c6515278c0    Liver Cancer - RIKEN, JP        2       1440.0  deceased       67.0    I       male    Liver-HCC
```

The script can also be run to include the actual mutant and wildtype transcripts in the table as well by using the ```-include-sequences``` flag:

```
➜ python seq-to-pheno/seq_to_pheno/tcga/create_protein_metadata_table.py -include-sequences
Metadata table saved to seq_to_pheno/tcga/data/protein_sequences_metadata.tsv
```

Inside table:
```
aliquot_id      mutated_protein_path    wildtype_protein_path   wgs_aliquot_id  Cancer Type     Cancer Stage    Donor Survival Time     Donor Vital Status      Donor Age at Diagnosis Tumour Grade    Donor Sex       Histology Abbreviation
8fb6f8bc-ddb8-11e4-ad8f-5ed8e2d07381    MLTRLQVLTLALFSKGFLLSLGDHNFLRREIKIEGDLVLGGLFPINEKGTGTEECGRINEDRGIQRLEAMLFAIDEINKDDYLLPGVKLSVHILDTCSRDTYALEQSLEFVRASLTKVDEAEYMCPDGSYAIQENIPLLIAGVIGGSYSSVSIQVANLLRLFQIPQISYASTSAKLSDKSRYDYFARTVPPDFYQAKAMAEILRFFNWTYVSTVASEGDYGETGIEAFEQEARLRNICIATAEKVGRSNIRKSYDSVIRELLQKPNARVVVLFMRSDDSRELIAAASRANASFTWVASDGWGAQESIIKGSEHVAYGAITLELASQPVRQFDRYFQSLNPYNNHRNPWFRDFWEQKFQCSLQNKRNHRRVCDKHLAIDSSNYEQESKIMFVVNAVYAMAHALHKMQRTLCPNTTKLCDAMKILDGKKLYKDYLLKINFTGADDNHVHLCQPEWLCGLGLFVCTQGSHHPVSTPEECCHTQTAPQQVQCQWNWDHILSVLCKHVCANGVQWAGSPRLHHLISVIVNCSSVLVFLDC*     MLTRLQVLTLALFSKGFLLSLGDHNFLRREIKIEGDLVLGGLFPINEKGTGTEECGRINEDRGIQRLEAMLFAIDEINKDDYLLPGVKLGVHILDTCSRDTYALEQSLEFVRASLTKVDEAEYMCPDGSYAIQENIPLLIAGVIGGSYSSVSIQVANLLRLFQIPQISYASTSAKLSDKSRYDYFARTVPPDFYQAKAMAEILRFFNWTYVSTVASEGDYGETGIEAFEQEARLRNICIATAEKVGRSNIRKSYDSVIRELLQKPNARVVVLFMRSDDSRELIAAASRANASFTWVASDGWGAQESIIKGSEHVAYGAITLELASQPVRQFDRYFQSLNPYNNHRNPWFRDFWEQKFQCSLQNKRNHRRVCDKHLAIDSSNYEQESKIMFVVNAVYAMAHALHKMQRTLCPNTTKLCDAMKILDGKKLYKDYLLKINFTGADDNHVHLCQPEWLCGLGLFVCTQGSHHPVSTPEECCHTQTAPQQVQCQWNWDHILSVLCKHVCANGVQWAGSPRLHHLISVIVNCSSVLVFLDC*     15fd8dc8-c622-11e3-bf01-24c6515278c0    Liver Cancer - RIKEN, JP       4       1200.0  deceased        62.0    III     male    Liver-HCC
```

## Next Steps

With this table we are at the point where we are very close to being ready to generate embeddings. We just need to create a way to subset the table intelligentl.

### Example Analysis:

- Load the seq_to_pheno/tcga/data/protein_sequences_metadata.tsv table into a dataframe
- Filter for a cancer type (Liver Cancer)
- Filter for a certain cancer stage (2 or T3N0MX)
- Order table by donor survival time
- Generate embeddings of wildtype and mutant transcripts.
- Experiment: 
  - Wildtype embeddings should cluster separately from mutant embeddings.
  - Mutant embeddings should cluster based on the extremes of donor survival time. 
  - Transcripts with the highest relational effect should show patterns across all/most subjects.

Don't feel committed to this exact example analysis. There are probably many ways we can create a reasonable and reproducable experiment. 