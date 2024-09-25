Here you will find a tutorial to perform a complete analysis for HTG Edge results and RNAseq with HTGAnalyzer.

Having the results of HTG Edge Reveal or RNAseq, this package will help you perform Quality Control, Differential Expression Analysis (DEA), Gene Set Enrichment Analysis (GSEA), Tumor Microenvironment Analysis (TME), and Survival Analysis.

# 1. INSTALATION:
To use this R package you have to perform the following code.
```{r}
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

library(devtools)
install_github("ISGLOBAL-Rakislova-Lab/HTGAnalyzer")
library (HTGAnalyzer)
```
# 2. TUTORIAL.
Once installed, you can start with the tutorial. This package, HTG_Analyzer, was originally designed for analyzing HTG Edge files. However, following the closure of the company, only certain quality control features were fully disclosed, particularly those not related to the transcriptomic panel. To address this, our package provides a simplified approach for users who may not have a background in bioinformatics. The primary function, HTG_auto, automates the entire analysis process—from data import and quality control to detailed analyses, complete with plots and tables. Every aspect is customizable, allowing users to easily modify and tailor the analyses to their specific needs. Although it was initially created for HTG, the package is also compatible with RNAseq data. Additionally, there are individual functions like HTG_DEA, HTG_survival, HTG_QC, and HTG_TME that offer more flexibility and precision for advanced users looking to perform specific analyses.

## 2.1. QUICK START.

HTGAnalyzer has a functions called HTG_auto which provide an easy way to perform all the analyses. Let's start by discussing HTG_auto.

### 2.1.1 HTG_auto.
The `HTG_auto` function in HTGAnalyzer automates a comprehensive analysis pipeline for HTG data. It integrates various analyses including Quality Control (QC), Differential Expression Analysis (DEA), Gene Set Enrichment Analysis (GSEA), Tumor Microenvironment Analysis (TME), and Survival Analysis. This function is designed to simplify the process, with default settings for ease of use, while also providing flexibility to modify parameters according to specific needs.
  
EXAMPLE:
imagine that your AnnotData looks like this: 

|       id       | HPV_status | Ciclina_D1 | FIGO_2021_STAGE | Recurrence | Recurrence_01 | Time_to_death_surv |
|----------------|------------|------------|-----------------|------------|----------------|---------------------|
| pacient_1      | Positive   | 10         | IIIB            | yes        | 1              | 1287                |
| pacient_2      | Positive   | 50         | II              | no         | 0              | 510                 |
| pacient_3      | Positive   | NA         | IIIC            | no         | 0              | 762                 |
| pacient_4      | Positive   | 0          | IB              | yes        | 1              | 1164                |
| pacient_5      | Negative   | 50         | IB              | yes        | 1              | 5844                |
| pacient_6      | Negative   | 80         | IB              | no         | 0              | 1436                |
| pacient_7      | Negative   | 20         | IB              | no         | 0              | 2145                |
| pacient_8      | Negative   | 50         | IB              | yes        | 1              | 2458                |
| pacient_9      | Negative   | NA         | IB              | yes        | 1              | 234                 |
| pacient_10     | Negative   | 70         | II              | no         | 0              | 996                 |
| pacient_11     | Negative   | 40         | IB              | no         | 0              | 768                 |

Then:
```{r}
FOR HTG:
HTG_auto <- function("~/counts.xlsx",
                     file_type = "HTG",
                     "~/AnnotData.xlsx",
                     design_formula = "HPV_status",
                     QC = TRUE,
                     heatmap_columns = c("HPV_status", "Ciclina_D1"),
                     contrast = c("HPV_status", "Positive", "Negative"),
                     variable_01 = "Recurrence_01",
                     time = "Time_to_death_surv")


FOR RNAseq:
HTG_auto <- function("~/counts.xlsx",
                     file_type = "RNAseq",
                     "~/AnnotData.xlsx",
                     design_formula = "HPV_status",
                     QC = TRUE,
                     heatmap_columns = c("HPV_status", "Ciclina_D1"),
                     contrast = c("HPV_status", "Positive", "Negative"),
                     variable_01 = "Recurrence_01",
                     time = "Time_to_death_surv")
```

For HTG data, setting QC = TRUE is essential to ensure quality control, which identifies and manages outliers before further analysis. With quality control enabled, the function performs Differential Expression Analysis (DEA) based on HPV_status, comparing positive and negative samples to identify significantly differentially expressed genes.

Following DEA, a heatmap will be generated using HPV_status and Ciclina_D1, providing insights into how well the data separates based on these variables. This visualization helps in understanding key distinguishing factors between the groups.

The DEA results will then be used for Gene Set Enrichment Analysis (GSEA) to determine if specific gene sets are enriched in differentially expressed genes, offering deeper biological insights.

Additionally, DEA findings will inform Tumor Microenvironment (TME) analysis, shedding light on tumor-environment interactions and their impact on gene expression and tumor biology.

Finally, a survival analysis will be conducted using the top genes from DEA and survival variables Recurrence_01 and Time_to_death_surv. This analysis evaluates the impact of gene expression on patient outcomes, providing valuable prognostic information.


## 2.2 IN-DEPTH GUIDE

For this part of the tutorial, we will combine principal and secondary functions to explore their capabilities. 

### 2.2.1 PRINCIPAL FUNCTIONS AND SECUNDARY FUNCTIONS 
Just in case you need more control of your data or you want to perform other analysis. HTGAnalyzer have principal and secundary fucntions to help reach your demands.

PRINCIPAL: 
* `HTG_import_counts`:Import counts data from an Excel file. It can be either HTG excel or RNAseq. The first row in the Excel file must contain the column headers, with "id" as the first column header followed by the names of each sample. Please note that you might need to modify the Excel file to ensure it is in the correct format for importing data into R. Ensure column names are free of special characters (e.g., spaces, (,), ?, `, ^, ., -, *, and others) to avoid import issues. The program will change the spaces into `_`
* `HTG_QC`: This function performs comprehensive quality control (QC) for HTG EdgeSeq data, including checks on positive values, library size, and various thresholds (e.g., negative control, genomic DNA, ERCC). It generates a data frame with summary statistics and sample ratios, highlights outliers with an optional heatmap, and saves plots in the working directory.
* `HTG_analysis`: This function conducts a comprehensive analysis pipeline including DESeq2 differential expression analysis (DEA), Gene Set Enrichment Analysis (GSEA), TME analysis, and survival analysis. The pipeline supports optional steps for generating volcano plots and heatmaps. The function is suitable for both HTG and RNA-seq data.

SECUNDARY:
* `HTG_subset`: This function subsets a data frame based on a specified prefix in the row names. It allows you to extract specific rows that match a given pattern and optionally normalizes the data using TPM (Transcripts Per Million) before subsetting. The function also provides the dimensions of the resulting data frame.
* `HTG_plotPCA`: This function generates a PCA plot specifically for genes, along with plots showing explained variability and accumulated explained variability. It highlights the samples that are the most distant from the center of the PCA plot.
* `quant_to_qual`: This function converts quantitative columns in a data frame into qualitative ones based on a specified threshold.

### 2.2.2 DATA IMPORT.

#### 2.2.2.1 HTG_import_counts

Let's start with the `HTG_import_counts` function. This function is designed to import count data (either RNAseq or HTG) into R from an Excel file. 

#### Typical HTG EdgeSeq Data Format

Results from the HTG EdgeSeq machine usually have a header that looks like this:

| Sample Name   | patient_1     | patient_2     | patient_3     | patient_4     | patient_5     | patient_6     | patient_7     | patient_8     |
|---------------|---------------|---------------|---------------|---------------|---------------|---------------|---------------|---------------|
| Sample ID     | 1             | 2             | 3             | 4             | 5             | 6             | 7             | 8             |
| Well          | A1            | B1            | C1            | D1            | E1            | F1            | G1            | H1            |
| Date Parsed   | 14/08/2024    | 14/08/2024    | 14/08/2024    | 14/08/2024    | 14/08/2024    | 14/08/2024    | 14/08/2024    | 14/08/2024    |
| Total Counts  | 9650628       | 33950226      | 24680689      | 10605963      | 7395173       | 4057094       | 20049502      | 13607678      |
| A3GALT2       | 0             | 0             | 164           | 19            | 46            | 1             | 0             | 0             |
| A4GALT        | 2             | 662           | 665           | 830           | 748           | 89            | 256           | 333           |
| A4GNT         | 0             | 1041          | 103           | 0             | 0             | 2             | 0             | 0             |
| AAAS          | 0             | 368           | 70            | 634           | 707           | 425           | 587           | 872           |
| AACS          | 1             | 0             | 562           | 728           | 250           | 209           | 365           | 445           |
| AADAC         | 0             | 382           | 0             | 5             | 144           | 0             | 0             | 1             |

The `HTG_import_counts` function can import this file, recognize the header, remove it, and keep the first row as column names and the first column as row names.

#### Alternative Data Format

If your Excel file has sample names as the first row and genes in the first column, like the format shown below, it will work as well:

| id      | patient_10 | patient_11 | patient_5 | patient_6 | patient_4 | patient_1 | patient_2 | patient_3 | patient_7 | patient_9 | patient_8 |
|---------|------------|------------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|
| A3GALT2 | 0          | 0          | 164       | 19        | 46        | 1         | 0         | 0         | 0         | 0         | 0         |
| A4GALT  | 2          | 662        | 665       | 830       | 748       | 89        | 256       | 333       | 210       | 2         | 504       |
| A4GNT   | 0          | 1041       | 103       | 0         | 0         | 2         | 0         | 0         | 84        | 0         | 0         |
| AAAS    | 0          | 368        | 70        | 634       | 707       | 425       | 587       | 872       | 542       | 0         | 370       |
| AACS    | 1          | 0          | 562       | 728       | 250       | 209       | 365       | 445       | 62        | 46        | 785       |
| AADAC   | 0          | 382        | 0         | 5         | 144       | 0         | 0         | 1         | 260       | 683       | 77        |

```{r}
counts<- HTG_import("path_to_HTG_database.xlsx")
head(counts)
```
After importing, your data will look like this:

|         | patient_10 | patient_11 | patient_5 | patient_6 | patient_4 | patient_1 | patient_2 | patient_3 | patient_7 | patient_9 | patient_8 |
|---------|------------|------------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|
| A3GALT2 | 0          | 0          | 164       | 19        | 46        | 1         | 0         | 0         | 0         | 0         | 0         |
| A4GALT  | 2          | 662        | 665       | 830       | 748       | 89        | 256       | 333       | 210       | 2         | 504       |
| A4GNT   | 0          | 1041       | 103       | 0         | 0         | 2         | 0         | 0         | 84        | 0         | 0         |
| AAAS    | 0          | 368        | 70        | 634       | 707       | 425       | 587       | 872       | 542       | 0         | 370       |
| AACS    | 1          | 0          | 562       | 728       | 250       | 209       | 365       | 445       | 62        | 46        | 785       |
| AADAC   | 0          | 382        | 0         | 5         | 144       | 0         | 0         | 1         | 260       | 683       | 77        |

This function can handle both formats efficiently, ensuring that your data is correctly imported into R for further analysis.

For AnnotData importation, you can easy do it by:
```{r}
AnnotData<- read_excel("path/to/your/annot_file.xlsx")
head(AnnotData)
```
We recomend you to avoir special characters on excel files (e.g., spaces, (,), ?, `, ^, ., -, *, and others)  to avoid analysis issues.

### 2.2.3 QUALITY CONTROL.

#### 2.2.3.1 HTG_QC
This function performs various quality control (QC) checks tailored for the HTG EdgeSeq transcriptomic panel, though thresholds can be adjusted as needed. QC checks include:

* **QC0**: % of positive values > 4%.
* **QC1**: Library size > 7e+06.
* **QC2**: Negative control threshold > 0.045.
* **QC3**: Genomic DNA threshold > 0.02.
* **QC4**: ERCC threshold > 0.025.

The function generates a data frame (which can be saved as a .csv file, with a preview shown in the R console) with:
* Sum of each probe for each sample (total genes, positive, negative, gdna, ercc).
* Ratios for each sample.
* Sample sizes.
* A PCA column indicating samples furthest from the center.
Additionally, a statistical .csv is generated with columns for Min, Max, Mean, Median, Mode, SD, Variance, Range, Q1, Q3, IQR, Skewness, Kurtosis, Missing, and CV.

The function includes a heatmap to highlight potential outliers, which will be saved in a vector. Additionally, it performs a PCA analysis that identifies and highlights the specified number of samples that are furthest from the center. The PCA results include plots showing the explained variability and the cumulative explained variability.

```{r}
outliers<- HTG_QC(counts_data_tutorial)
outliers
```

#### 2.2.3.2 HTG_plotPCA and HTG_calculate_summary_stats
This QC process is designed primarily for HTG data, assuming that RNAseq data has already passed the necessary controls. However, should additional checks be needed, we use two secondary functions:

1. **`HTG_plotPCA`**: Performs Principal Component Analysis (PCA).
2. **`HTG_calculate_summary_stats`**: Calculates detailed summary statistics.

The summary statistics include columns for:
- **Min**: Minimum value
- **Max**: Maximum value
- **Mean**: Mean
- **Median**: Median
- **Mode**: Mode
- **SD**: Standard Deviation
- **Variance**: Variance
- **Range**: Range (Max - Min)
- **Q1**: First Quartile
- **Q3**: Third Quartile
- **IQR**: Interquartile Range (Q3 - Q1)
- **Skewness**: Skewness
- **Kurtosis**: Kurtosis
- **Missing**: Count of missing values
- **CV**: Coefficient of Variation (SD / Mean)

```{r}
# PCA FOR HTG
PCA_HTG<- HTG_plotPCA(counts_data_tutorial, n_samples = 0,pattern = "^NC-|^POS-|^GDNA-|^ERCC-")
PCA_HTG
summary <- HTG_calculate_summary_stats(counts_data_tutorial, pattern = "^NC-|^POS-|^GDNA-|^ERCC-")
summary


# PCA FOR RNAseq
PCA_RNAseq<- HTG_plotPCA(counts_data_tutorial, n_samples = 0,pattern = NULL)
PCA_RNAseq
summary <- HTG_calculate_summary_stats(counts_data_tutorial)
summary
```

#### 2.2.3.3 HTG_subset
The `HTG_subset` function is a versatile tool for quickly extracting specific genes or probes from your data, whether it be raw counts or results from differential expression analysis. This functionality is especially useful for focusing on individual genes of interest and facilitating detailed visualization.

EXAMPLES:

```r
# Subset ERCC probes from counts data
  ERCC <- HTG_subset(counts_data_tutorial, "ERCC")

# Extracting a Specific Gene (AAAS) from Counts Data
# Subset AAAS gene from counts data without normalization
counts_AAAS <- HTG_subset(counts_data_tutorial, "AAAS")

# Subset AAAS gene from counts data with normalization
counts_AAAS_normalize <- HTG_subset(counts_data_tutorial, "AAAS", normalize = TRUE)

# Extracting a Specific Gene (AAAS) from Differential Expression Results
# Subset AAAS gene from differential expression results
res_AAAS <- HTG_subset(res_tutorial, "AAAS")
```

### 2.2.4 ANALYSIS
After completing the quality control (QC) steps, you can perform up to four different types of analysis:

* Differential Expression Analysis
* Gene Set Enrichment Analysis (GSEA)
* Tumor Microenvironment Analysis
* Survival Analysis

Since HTG data may have outliers while RNA-seq data typically does not, we will provide two examples for each type of analysis specific to HTG and RNA-seq datasets. Please note that the documentation will include various adjustable parameters for each function, allowing you to tailor the analysis to your specific needs.


#### 2.2.4.5 HTG_analysis
All these analyses can be performed using the `HTG_analysis` function.  This function facilitates the execution of differential expression analysis, Gene Set Enrichment Analysis (GSEA), tumor microenvironment analysis, and survival analysis. In the two examples provided in this tutorial, the function is configured to perform all analyses. However, you can customize the function to execute only the analyses you require by setting the relevant parameters to `TRUE` or `FALSE`.


```{r}
# EXAMPLE HTG:
ALL_analysis <- HTG_analysis(
  outliers = outliers, 
  pattern = "^NC-|^POS-|^GDNA-|^ERCC-", 
  counts_data = counts_data_tutorial, 
  col_data = AnnotData_tutorial, 
  design_formula = "HPV_status",
  heatmap_columns = c("HPV_status", "Ciclina_D1"), 
  contrast = c("HPV_status", "Positive", "Negative"), 
  pCutoff = 5e-2, 
  variable_01 = "Recurrence_01", 
  time = "Time_to_death_surv", 
  GSEA = TRUE, 
  survival_analysis = TRUE)


# EXAMPLE RNAseq:
ALL_analysis <- HTG_analysis(
  counts_data = counts_data_tutorial, 
  col_data = AnnotData_tutorial, 
  design_formula = "HPV_status",
  heatmap_columns = c("HPV_status", "Ciclina_D1"), 
  contrast = c("HPV_status", "Positive", "Negative"), 
  pCutoff = 5e-2, 
  variable_01 = "Recurrence_01", 
  time = "Time_to_death_surv", 
  GSEA = TRUE, 
  survival_analysis = TRUE)
```
The `HTG_analysis` function generates several Excel files and plots. These can be visualized and saved as PDFs if needed.


#### 2.2.4.2 HTG_DEA

Each analysis can be performed separately, giving you greater control and flexibility over the process. For example, the `HTG_DEA` function allows you to tailor the differential expression analysis to meet specific needs and datasets. It provides options to apply filters and perform log-fold change (LFC) shrinkage as required, enabling you to customize the analysis according to your objectives.


```{r}
# EXAMPLE FOR HTG: 
results <- HTG_DEA(
  outliers = outliers,
  counts_data = counts_data_tutorial,
  col_data = AnnotData_tutorial,
  design_formula = "HPV_status",
  heatmap_columns = c("HPV_status", "Ciclina"),
  contrast = c("HPV_status", "Positive", "Negative"),
  pattern = "^NC-|^POS-|^GDNA-|^ERCC-",
  percentage_gene = 0.2,
  percentage_zero = 0.2,
  threshold_gene = 200,
  threshold_subject = 10,
  pCutoff = 5e-2,
  apply_filtering = TRUE,
  apply_lfc_shrinkage = TRUE,
  extract_shrinkage = TRUE
)

# EXAMPLE FOR RNAseq:
results <- HTG_DEA(
  counts_data = counts_data_tutorial,
  col_data = AnnotData_tutorial,
  design_formula = "HPV_status",
  heatmap_columns = c("HPV_status", "Ciclina"),
  contrast = c("HPV_status", "Positive", "Negative"),
  pattern = "^NC-|^POS-|^GDNA-|^ERCC-",
  remove_outliers = FALSE,
  percentage_gene = 0.2,
  percentage_zero = 0.2,
  threshold_gene = 200,
  threshold_subject = 10,
  pCutoff = 5e-2,
  apply_filtering = TRUE,
  apply_lfc_shrinkage = TRUE,
  extract_shrinkage = TRUE
)
```
#### 2.2.4.3 HTG_GSEA

If you have already conducted differential expression analysis using other packages, you can still perform Gene Set Enrichment Analysis (GSEA) with the `HTG_GSEA` function.

```{r}
HTG_GSEA(res_tutorial)
```
#### 2.2.4.4 HTG_TME

The `HTG_TME` function is designed to perform tumor microenvironment (TME) analysis using three different methods. Unlike the HTG_analysis function, which provides results in various formats, `HTG_TME` returns an object containing three separate tables for each method. This function is specifically suited for use with differential expression data from the DESeq2 library.

```{r}
# FOR HTG
TME_data<- HTG_TME(outliers = outliers_tutorial, pattern =  "^NC-|^POS-|^GDNA-|^ERCC-", counts_data = counts_data_tutorial , AnnotData= AnnotData_tutorial, design_formula = "HPV_status")
TME_data$EPIC
TME_data$QTI
TME_data$XCELL


# FOR RNAseq
TME_data<- HTG_TME(counts_data = counts_data_tutorial, AnnotData = AnnotData_tutorial, design_formula = "HPV_status")
TME_data$EPIC
TME_data$QTI
TME_data$XCELL
```

#### 2.2.4.5 HTG_survival

The `HTG_survival` function is designed for performing survival analysis and offers flexibility for various use cases:

* Top 10 Genes from Differential Expression Analysis: Conduct survival analysis on the top 10 genes identified from your differential expression results.
* Custom Genes: Perform survival analysis on a selected set of genes that you specify.
* Tumor Microenvironment Analysis Results: Analyze survival data based on results from tumor microenvironment studies.


```{r}

# Survival from Top10 genes
survival_res<- HTG_survival(variable_01 = "Recurrence_01", time = "Time_to_death_surv", col_data = AnnotData_tutorial,
              count_data = counts_data_tutorial, res = res_tutorial, genes_to_use = NULL,
              outliers = outliers_tutorial, pattern = "^NC-|^POS-|^GDNA-|^ERCC-", remove_outliers = TRUE)

# Survial of CCND1 gene
survival_gene_to_use<- HTG_survival(variable_01 = "Recurrence_01", time = "Time_to_death_surv", col_data = AnnotData_tutorial,
              count_data = counts_data_tutorial, res = res_tutorial, genes_to_use = c("LCP1", "OMA1"),
              outliers = outliers_tutorial, pattern = "^NC-|^POS-|^GDNA-|^ERCC-", remove_outliers = TRUE)

# Survival of EPIC TME results.
survival_TME<- HTG_survival(variable_01 = "Recurrence_01", time = "Time_to_death_surv", col_data = AnnotData_tutorial,
              count_data = counts_data_tutorial, res = NULL, genes_to_use = NULL, TME = TME_data_tutorial$EPIC,
              outliers = outliers_tutorial, pattern = "^NC-|^POS-|^GDNA-|^ERCC-", remove_outliers = TRUE)


FOR RNAseq: 
# Survival from Top10 genes
survival_res<- HTG_survival(variable_01 = "Recurrence_01", time = "Time_to_death_surv", col_data = AnnotData_tutorial, count_data = counts_data_tutorial, res = res_tutorial)
                         
# Survial of CCND1 gene
survival_gene_to_use<- HTG_survival(variable_01 = "Recurrence_01", time = "Time_to_death_surv", col_data = AnnotData_tutorial, count_data = counts_data_tutorial, genes_to_use = "CCND1")

# Survival of EPIC TME results.
survival_TME<- HTG_survival(variable_01 = "Recurrence_01", time = "Time_to_death_surv", col_data = AnnotData_tutorial, count_data = counts_data_tutorial, TME = TME_data_tutorial$EPIC)
```
