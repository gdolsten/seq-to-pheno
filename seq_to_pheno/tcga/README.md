# TCGA

Space for processing TCGA data

## Datasets
The ICGC (International Cancer Genome Consortium) in collaboration with TCGA and the Sanger Institute has created a hub for a wide-range of NGS data and metadata across various cancer types.

Instructions on how to access the data are [found on this page](https://docs.icgc-argo.org/docs/data-access/icgc-25k-data).

### Get
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
➜ head -5 rnaseq.extended.metadata.aliquot_id.V4.tsv | cut -f1,3,76,48 
aliquot_id      project_code    donor_survival_time     GTex.tissue
b337121c-9821-4644-820e-b8c477f6c38a    GBM-US  489     Brain
612ef912-5a28-4c11-8703-3376f51afef5    GBM-US  317     Brain
56a705f4-fd28-44ff-8a3c-85bc4300c760    GBM-US  737     Brain
28239a0e-9990-49ef-a159-ba63fb078c77    GBM-US  485     Brain
```
