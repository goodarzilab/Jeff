# BRCA oncRNA Project
In this project, we identify breast cancer specific orphan non-coding RNAs (oncRNAs) in our cancer cell lines and TCGA BRCA smRNA datasets. We also include downstream analysis to prioritize oncRNAs for experimental validation.
All scripts run can be found in the `scripts/` directory.

## 1) Preprocessing
All preprocessing steps of the raw data can be found in the following notebook:
- `preprocess_data.ipynb`

Datasets used in this project include:
- small RNAseq data of in-house breast cancer and normal breast tissue cell lines.
- TCGA BRCA (breast invasice carcinoma) smRNA datasets and all TCGA normal tissue smRNA datasets.
- exRNA Atlas non-cancerous exosomal, small RNAseq datasets.

## 2) Cancer-Enriched Loci Discovery
For this discovery step, we implemented three different methods analysis to identify potential orphan non-coding RNAs. \
These approaches include: Fisher Exact Test to identify cancer-specifc RNAs, DESeq to identifiy cancer-enriched RNAs, and individual DESeq 
analysis based on individual subtypes. <br>
Code can be found in the following notebooks:
- `cancer_enriched_loci_fisher_analysis.ipynb`
- `cancer_enriched_loci_DE_subtype_analysis.ipynb`
- `cancer_enriched_loci_DE_analysis.ipynb`

## 3) BRCA Normal Thresholding
Here we compare our breast cancer cell lines enriched RNAs with TCGA BRCA samples and and all TCGA normal samples. We perform a thresholding experiment to set a filter threshold (for presence in TCGA normal samples) for our putative oncRNAs. This is done to narrow down our pool of putative oncRNAs under consideration by eliminating RNAs found in abundance in normal tissues before statistical testing (Fisher Exact Test) for significance in TCGA BRCA. Resulting oncRNAs that are found to be significantly enriched in our cell lines and in TCGA BRCA are plotted with a heatmap. <br>
Code can be found in the notebook:
- `BRCA_normal_thresholding.ipynb`

## 4) Survival Analysis and Metastatic DE Analysis
- `survival_analysis.ipynb`
- `metastatic_cell_lines_DE.ipynb`

## 5) Finalizing oncRNAs for Experimental Validation
Here we prioritize oncRNAs identified in the Nature Med cell lines and pan-cancer BRCA oncRNAs observed in the miniENCODE cell lines to move into experimental validation. 
- `oncRNAs_combined.ipynb` 


Here we trim and finalize the oncRNA loci and generate the sequences for testing.
- `oncRNAs_final_list.ipynb`

