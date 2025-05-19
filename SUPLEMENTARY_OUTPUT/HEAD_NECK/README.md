#TCGA Head and Neck Squamous Cell Carcinoma (HNSC) – Gene Expression and Survival Analysis

This repository contains results from the analysis of gene expression and survival in the **TCGA Head and Neck Squamous Cell Carcinoma (TCGA-HNSC)** cohort using HTG_auto function from HTGAnalyzer package.  
Data were accessed through [cBioPortal](https://www.cbioportal.org/study/summary?id=hnsc_tcga_gdc), originally released by the Genomic Data Commons (GDC) in **July 2024**.

##Dataset Overview

- **Samples**: 502 primary tumor samples  
- **Sex distribution**: 133 female, 369 male patients  
- **Included data**:
  - Gene expression (raw read counts)
  - Harmonized clinical annotations
- **Key clinical variables**:
  - `SEX` used as a contrast factor in analyses
  - `OS_months` and `OS_status` used to define overall survival

Gene identifiers were standardized by converting Entrez IDs to HGNC symbols using the `org.Hs.eg.db` database to ensure consistency across datasets.

---

## Files

### 📄 `HEAD_AND_NECK_RESULTS.pdf`
Summary report with key findings from the analysis.

### 📄 `imm_epic.csv`
Immune infiltration estimates using the **EPIC** method.

### 📄 `imm_qti.csv`
Immune composition estimates using the **QTi** approach.

### 📄 `imm_xcell.csv`
Tumor microenvironment estimates based on the **xCell** method.

---

### 📄 Differential Survival Analysis per Gene

These files contain results from survival and differential expression analysis stratified by sex, focused on Y-chromosome genes:

- `surv_diff_summary_EIF1AY_mRNA_expression.csv`
- `surv_diff_summary_KDM5D_mRNA_expression.csv`
- `surv_diff_summary_RPS4Y1_mRNA_expression.csv`
- `surv_diff_summary_TXLNGY_mRNA_expression.csv`
- `surv_diff_summary_USP9Y_mRNA_expression.csv`
- `surv_diff_summary_UTY_mRNA_expression.csv`
- `surv_diff_summary_ZFY_mRNA_expression.csv`

Each file summarizes mRNA expression values and survival differences between male and female patients for a specific gene.

---

## Citation

If you use this repository or data in your work, please cite the original sources:

> TCGA Research Network: https://www.cancer.gov/tcga  
> cBioPortal for Cancer Genomics: https://www.cbioportal.org/

---

## Contact

For questions or collaborations
