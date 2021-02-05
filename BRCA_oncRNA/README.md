# BRCA oncRNA Project
In this project, we identify breast cancer specific orphan non-coding RNAs (oncRNAs) in our cancer cell lines and TCGA BRCA smRNA datasets. We also include downstream analysis to prioritize oncRNAs for experimental validation.
All scripts run can be found in the `scripts/` directory.

## 1) Preprocessing
All preprocessing steps of the raw data can be found in the following notebook:
- `preprocess_data.ipynb`

Datasets used in this project include:
- small RNAseq data of in-house breast cancer and normal breast tissue cell lines.
- TCGA BRCA smRNA datasets and all TCGA normal tissue smRNA datasets.
- exRNA Atlas non-cancerous exosomal, small RNAseq datasets.

## 2) Cancer-Enriched Loci Discovery
For this discovery step, we implemented three different methods analysis to identify potential orphan non-coding RNAs. \
These approaches include: Fisher Exact Test to identify cancer-specifc RNAs, DESeq to identifiy cancer-enriched RNAs, and individual DESeq 
analysis based on individual subtypes. <br>
Code can be found in the following notebook
- `cancer_enriched_loci_fisher_analysis.ipynb`
- `cancer_enriched_loci_DE_subtype_analysis.ipynb`
- `cancer_enriched_loci_DE_analysis.ipynb`

## 3) BRCA Normal Thresholding


## 4) Survival Analysis


## 5) RNA Annotations







