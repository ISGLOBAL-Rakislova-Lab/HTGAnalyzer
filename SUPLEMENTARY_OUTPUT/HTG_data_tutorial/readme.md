## Files

This project analyzes transcriptomic data generated with the **HTG Transcriptome Panel**, including **Differential Expression Analysis (DEA)**, **Gene Set Enrichment Analysis (GSEA)**, **Tumor Microenvironment Estimation (TME)**, and **Survival Analysis**. The analysis was performed using the `HTG_auto()` pipeline.

### 🧪 HTG Transcriptome Panel

A total of 58 FFPE vulvar squamous cell carcinoma (VSCC) samples were processed using the HTG Transcriptome Panel, which includes:

- 19,398 probes for human RNA transcripts  
- 218 control probes for QC  
  - 4 positive controls  
  - 100 negative controls  
  - 22 genomic DNA probes  
  - 92 ERCC probes  

Samples were processed using the **HTG EdgeSeq Processor** and sequenced on an **Illumina platform**. QC was assessed on 50 of the 58 samples using the **HTG EdgeSeq Reveal** software. See Supplementary Table S1 for detailed QC results.

---

### 🔧 Code Execution

The full pipeline was executed using the `HTG_auto()` function in R. The exact command can be found in the project's scripts.
```{r}
 HTG_auto(
    "~/counts.xlsx",
    file_type        = "HTG",
    "~/AnnotData.xlsx",
    design_formula   = "HPV_status",
    QC               = TRUE,
    heatmap_columns  = c("HPV_status","Cyclin_D1"),
    contrast         = c("HPV_status","Associated","Independent"),
    variable_01      = "Recurrence_01",
    time             = "time_to_recurrence",
    DEA              = TRUE,
    genes_to_use     = NULL,
    remove_outliers  = TRUE,
    GSEA             = TRUE,
    generate_heatmap = TRUE,
    TME              = TRUE,
    survival_analysis= TRUE
  )
```

---

### 📁 Raw Data
You can download the **mock raw data** used in the tutorial here:  
🔗 [HTG_data_tutorial](https://github.com/ISGLOBAL-Rakislova-Lab/HTGAnalyzer/tree/main/SUPLEMENTARY_OUTPUT/HTG_data_tutorial/HTG_data_tutorial)

Contents:
- **`count.xlsx`** – Raw gene expression matrix (HTG format, integers only)  
- **`AnnotData.xlsx`** – Sample annotation / clinical data  

---

### 🧬 Differential Expression Analysis (DEA)

- `results_DEA.csv`: Differentially expressed genes  
- `DEA_plots_HTG_analysis.pdf`: Plots of top DEGs and volcano plots

---

### 📈 Gene Set Enrichment Analysis (GSEA)

- `gene_list.csv`: Ranked gene list used for GSEA  
- `go_results.csv`: GO terms significantly enriched  
- `kegg_gene_list.csv`: Ranked list for KEGG analysis  
- `GSEA_analysis_plots1_of_2.pdf` & `GSEA_analysis_plots2_of_2.pdf`: GSEA visualizations  

---

### 🔥 Heatmaps

- `HEATMAP_analysis.pdf`: Heatmap of top DE genes  
- `HTG_QC_heatmap.pdf`: QC metrics heatmap

---

### 📊 Quality Control (QC)

- `QC_results.csv`: QC summary of all samples  
- `QC_plots.pdf`: QC metric distributions  
- `QC_plots_violin_plot.pdf`: Violin plots of QC parameters  

---

### 🌿 Tumor Microenvironment (TME)

- `imm_epic.csv`, `imm_qti.csv`, `imm_xcell.csv`: Cell-type deconvolution matrices  
- `plot_cell_fraction_Average_cell_fraction_EPIC.pdf`  
- `plot_cell_fraction_Average_cell_fraction_quantTiseq.pdf`  
- `plot_cell_fraction_Average_cell_fraction_xCell.pdf`  
- `plots_TME_HeatmapEPIC.pdf`, `plots_TME_Heatmap_qti.pdf`, `plots_TME_Heatmap_xcell.pdf`  
- `plots_imm_EPIC.pdf`, `plots_imm_qti.pdf`, `plots_imm_xcell.pdf`: Barplots of cell-type proportions  

---

### ⏳ Survival Analysis

- `summary_stats.csv`: Clinical and survival summary of the cohort

---

## Citation

> Please cite the HTG EdgeSeq platform, and refer to the original dataset as described in the manuscript. For software and pipeline use, credit the authors of this repository and relevant R packages.

---

## Contact

For questions or contributions, please open an issue or contact the project maintainers.
