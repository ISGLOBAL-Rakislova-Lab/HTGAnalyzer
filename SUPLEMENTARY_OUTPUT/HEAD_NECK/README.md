## Files

This project followed a comprehensive analysis pipeline including **Differential Expression Analysis (DEA)**, **Gene Set Enrichment Analysis (GSEA)**, **Tumor Microenvironment Estimation (TME)**, and **Survival Analysis**. All steps were run using the `HTG_auto()` function.

### 🔧 Code

The pipeline was executed with the following command in R:

```r
HTG_auto(
  "~/counts_filtered.xlsx",
  file_type                 = "RNAseq",
  "~/head_neck_filtered.xlsx",
  pattern                   = NULL,
  QC                        = FALSE,
  threshold_superior_pos    = 5,
  threshold_line_pos        = 4,
  threshold_inferior_lib    = 5e+06,
  threshold_lib             = 7e+06,
  threshold_superior_nc     = 0.05,
  threshold_line_nc         = 0.045,
  threshold_superior_gdna   = 0.025,
  threshold_line_gdna       = 0.02,
  threshold_superior_ercc   = 0.03,
  threshold_line_ercc       = 0.025,
  threshold_inferior_median = 3,
  threshold_line_median     = 5,
  save_csv                  = FALSE,
  csv_file                  = "QC_results_2.csv",
  design_formula            = "SEX",
  percentage_gene           = 0.2,
  threshold_gene            = 200,
  threshold_subject         = 10,
  genes_to_use              = c("CDKN2A", "TP63", "IRF6", "KRT4"),
  heatmap_columns           = c("SEX", "ALCOHOL_HISTORY_DOCUMENTED"),
  contrast                  = c("SEX", "Male", "Female"),
  pCutoff                   = 0.05,
  variable_01               = "OS_STATUS",
  time                      = "OS_MONTHS",
  DEA                       = TRUE,
  remove_outliers           = FALSE,
  GSEA                      = TRUE,
  generate_heatmap          = TRUE,
  TME                       = TRUE,
  survival_analysis         = TRUE
)
```

---

### 📁 Raw Data

- `head_neck_filtered.xlsx`: Clinical metadata used in the analysis.
- `counts_filtered.xlsx`: Filtered gene expression matrix (raw read counts).

---

### 🧬 Differential Expression Analysis (DEA)

- `results_DEA.csv`: Differentially expressed genes by sex (`SEX` variable: Male vs Female).

---

### 📈 Gene Set Enrichment Analysis (GSEA)

- `gene_list.csv`: Ranked gene list used for GSEA, derived from DEA results.
- `HEAD_AND_NECK_RESULTS.pdf`: Summary report including visualizations and interpretations of GSEA results.

---

### 🌿 Tumor Microenvironment (TME) Analysis

- `imm_epic.csv`: Cell-type composition inferred using the **EPIC** method.
- `imm_qti.csv`: Immune cell abundance estimated by the **QTi** method.
- `imm_xcell.csv`: Microenvironmental profiling based on **xCell** scores.

---

### ⏳ Survival Analysis

These files contain survival and expression analyses for selected Y-chromosome genes, highlighting sex-based differences:

- `surv_diff_summary_EIF1AY_mRNA_expression.csv`
- `surv_diff_summary_KDM5D_mRNA_expression.csv`
- `surv_diff_summary_RPS4Y1_mRNA_expression.csv`
- `surv_diff_summary_TXLNGY_mRNA_expression.csv`
- `surv_diff_summary_USP9Y_mRNA_expression.csv`
- `surv_diff_summary_UTY_mRNA_expression.csv`
- `surv_diff_summary_ZFY_mRNA_expression.csv`

Each file includes survival metrics and expression levels stratified by sex.

---

## Citation

If you use this repository or data in your work, please cite the original sources:

> TCGA Research Network: https://www.cancer.gov/tcga  
> cBioPortal for Cancer Genomics: https://www.cbioportal.org/

---

## Contact

For questions or collaboration opportunities, feel free to open an issue or reach out to the repository maintainer.

