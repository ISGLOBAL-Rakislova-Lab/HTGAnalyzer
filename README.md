<!-- badges: start -->
[![R-CMD-check](https://github.com/ISGLOBAL-Rakislova-Lab/HTGAnalyzer/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ISGLOBAL-Rakislova-Lab/HTGAnalyzer/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->
[![Codecov test coverage](https://codecov.io/gh/ISGLOBAL-Rakislova-Lab/HTGAnalyzer/branch/master/graph/badge.svg)](https://app.codecov.io/gh/ISGLOBAL-Rakislova-Lab/HTGAnalyzer?branch=master)
[![License](https://img.shields.io/badge/license-GPL2-blue.svg)](https://github.com/ISGLOBAL-Rakislova-Lab/HTGAnalyzer/blob/master/LICENSE)
[![Documentation](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://ISGLOBAL-Rakislova-Lab.github.io/HTGAnalyzer)


# TUTORIAL

Here you will find a tutorial to perform a complete analysis for HTG Edge results and RNAseq with HTGAnalyzer.

Having the results of **HTG Edge Reveal** or **RNAseq**, this package will help you perform Quality Control, Differential Expression Analysis (DEA), Gene Set Enrichment Analysis (GSEA), Tumor Microenvironment Analysis (TME), and Survival Analysis.
![FIGURE_1](https://github.com/user-attachments/assets/49a47a9b-c6ae-4231-938f-9654c44995a3)

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
 **NOTE**: To install the `HTGAnalyzer` package, please ensure that the following GitHub packages are installed:

- [omnideconv/immunedeconv](https://github.com/omnideconv/immunedeconv)
- [dviraran/xCell](https://github.com/dviraran/xCell)
- [GfellerLab/EPIC](https://github.com/GfellerLab/EPIC)
- [IOBR/IOBR](https://github.com/IOBR/IOBR)
- [kevinblighe/EnhancedVolcano](https://github.com/kevinblighe/EnhancedVolcano)

Visit their GitHub pages for installation instructions.

# 2. TUTORIAL.
The HTGAnalyzer package is designed for users with limited bioinformatics experience who need a straightforward tool for transcriptomic data analysis. Initially created for HTG Edge files, the package now supports RNAseq data and addresses the gap left by the closure of HTG Edge. HTGAnalyzer provides an easy-to-use solution for conducting Quality Control (**QC**), Differential Expression Analysis (**DEA**), Gene Set Enrichment Analysis (**GSEA**), Tumor Microenvironment Analysis (**TME**), and **Survival Analysis** for transcriptomic panels.

## 2.0 DATA IMPUT:
To use this package, you'll need:
* An Excel file containing data counts, either from HTG or RNAseq,
* An Excel file with annotation or clinical data.

## 2.1. QUICK START.
HTGAnalyzer has a functions called `HTG_auto` which provide an easy way to perform all the analysis.

### 2.1.1 HTG_auto.
#### 2.1.1.1 HTG_auto: all the analysis
The `HTG_auto` function in HTGAnalyzer automates a comprehensive analysis pipeline for HTG data. It integrates various analysis including Quality Control (**QC**), Differential Expression Analysis (**DEA**), Gene Set Enrichment Analysis (**GSEA**), Tumor Microenvironment Analysis (**TME**), and **Survival Analysis**. This function is designed to simplify the process, with default settings for ease of use, while also providing flexibility to modify parameters according to specific needs.
  
**EXAMPLE:**
imagine that your AnnotData looks like this: 

** **IS VERY IMPORTANT THAT THE SAMPLE NAME IN AnnotData.xlsx HAS TO BE "id"** **

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
##### Replace "~/counts.xlsx" and "~/AnnotData.xlsx" with your authentic file paths.
##### You can find the example Excel files in the SUPPLEMENTARY_OUTPUT folder on GitHub.

# FOR HTG:
HTG_auto_HTG <- HTG_auto("~/counts.xlsx",
                     file_type = "HTG",
                     "~/AnnotData.xlsx",
                     design_formula = "HPV_status",
                     QC = TRUE,
                     heatmap_columns = c("HPV_status", "Ciclina_D1"),
                     contrast = c("HPV_status", "Positive", "Negative"),
                     variable_01 = "Recurrence_01",
                     time = "Time_to_death_surv",
                     DEA = TRUE,
                     remove_outliers = TRUE,
                     GSEA = TRUE,
                     generate_heatmap = TRUE,
                     TME = TRUE,
                     survival_analysis = TRUE)

# FOR RNAseq:
HTG_auto_RNA <- HTG_auto("~/counts.xlsx",
                     file_type = "RNAseq",
                     "~/AnnotData.xlsx",
                     design_formula = "HPV_status",
                     QC = FALSE, #it uses probes which are not in RNAseq data. 
                     heatmap_columns = c("HPV_status", "Ciclina_D1"),
                     contrast = c("HPV_status", "Positive", "Negative"),
                     variable_01 = "Recurrence_01",
                     time = "Time_to_death_surv",
                     DEA = TRUE,
                     remove_outliers = FALSE,
                     GSEA = TRUE,
                     generate_heatmap = TRUE,
                     TME = TRUE,
                     survival_analysis = TRUE)
```

The code generates a PDF and CSV with detailed analysis results. It performs Differential Expression Analysis (DEA) comparing HPV_status (positive vs. negative), removes outliers, and conducts Gene Set Enrichment Analysis (GSEA). A heatmap visualizes the separation based on HPV status and Ciclina_D1, with the option to add more columns. Additionally, it assesses the Tumor Microenvironment (TME) and performs survival analysis using top DEA genes, along with Recurrence_01 and Time_to_death_surv, to evaluate their impact on patient outcomes.

#### 2.1.1.2 HTG_auto:skiping analysis.
If you want to skip an analysis, simply set the value to FALSE.

```{r}

##### Replace "~/counts.xlsx" and "~/AnnotData.xlsx" with your authentic file paths.
##### You can find the example Excel files in the SUPPLEMENTARY_OUTPUT folder on GitHub.

# WITHOUT TME
HTG_auto_HTG <- HTG_auto("~/counts.xlsx",
                     file_type = "HTG",
                     "~/AnnotData.xlsx",
                     design_formula = "HPV_status",
                     QC = TRUE,
                     heatmap_columns = c("HPV_status", "Ciclina_D1"),
                     contrast = c("HPV_status", "Positive", "Negative"),
                     variable_01 = NULL,  #As you will not need this variable. You can keep it as NULL
                     time = "Time_to_death_surv",
                     DEA = TRUE,
                     remove_outliers = TRUE,
                     GSEA = TRUE,
                     generate_heatmap = TRUE,
                     TME = FALSE, # Keep it false. 
                     survival_analysis = TRUE)

# WITHOUT SURVIVAL
HTG_auto_HTG <- HTG_auto("~/counts.xlsx",
                     file_type = "HTG",
                     "~/AnnotData.xlsx",
                     design_formula = "HPV_status",
                     QC = FALSE,
                     heatmap_columns = c("HPV_status", "Ciclina_D1"),
                     contrast = c("HPV_status", "Positive", "Negative"),
                     variable_01 = "Recurrence_01",  
                     time = NULL, #As you will not need this variable. You can keep it as NULL
                     DEA = TRUE,
                     remove_outliers = TRUE,
                     GSEA = TRUE,
                     generate_heatmap = TRUE,
                     TME = TRUE,
                     survival_analysis = FALSE) # Keep it false
```

## 2.2 IN-DEPTH GUIDE
If you prefer to have more control over the analysis process or if you have already completed some analysis, and `HTG_auto` does not fully meet your needs, we offer additional options to tailor your workflow.

### 2.2.2 DATA IMPORT.
**NOTE**: We recomend you to avoir special characters on excel files (e.g., spaces, (,), ?, `, ^, ., -, *, and others)  to avoid analysis issues.

#### 2.2.2.1 HTG_import_counts
Let's start with the `HTG_import_counts` function. This function is designed to import count data (either **RNAseq or HTG**) into R from an Excel file. 

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

| id           | patient_1     | patient_2     | patient_3     | patient_4      | patient_5    | patient_6     | patient_7    | patient_8      |
|--------------|---------------|---------------|---------------|----------------|--------------|---------------|--------------|----------------|
| A3GALT2      | 0             | 0             | 164           | 19             | 46           | 1             | 0            | 0              |
| A4GALT       | 2             | 662           | 665           | 830            | 748          | 89            | 256          | 333            |
| A4GNT        | 0             | 1041          | 103           | 0              | 0            | 2             | 0            | 0              |
| AAAS         | 0             | 368           | 70            | 634            | 707          | 425           | 587          | 872            |
| AACS         | 1             | 0             | 562           | 728            | 250          | 209           | 365          | 445            |
| AADAC        | 0             | 382           | 0             | 5              | 144          | 0             | 0            | 1              |

```{r}

##### Replace "~/counts.xlsx" and "~/AnnotData.xlsx" with your authentic file paths.
##### You can find the example Excel files in the SUPPLEMENTARY_OUTPUT folder on GitHub.

# FOR HTG
counts<- HTG_import_counts("~/counts.xlsx", "HTG")
head(counts)

# FOR RNAseq
rna_data <- HTG_import_counts("~/counts.xlsx", "RNAseq")
head(rna_data)
```
After importing, your data will be ready for future steps. 

#### 2.2.2.2 AnnotData importation:
You can easy do it by:
```{r}
AnnotData<- read_excel("path/to/your/annot_file.xlsx")
head(AnnotData)
```

### 2.2.3 QUALITY CONTROL.
#### 2.2.3.1 HTG_QC
This function performs various quality control (QC) checks tailored for the HTG EdgeSeq transcriptomic panel, though thresholds can be adjusted as needed.
From now on, we will use the tutorial dataset to facilitate the execution.
```{r}
outliers<- HTG_QC(counts_data_tutorial)
outliers
```

#### 2.2.3.2 HTG_calculate_summary_stats
The `HTG_QC` function has already integrate this function. However, should additional checks be needed (for example for **RNAseq data** ,you can use `HTG_calculate_summary_stats` to calculates detailed summary statistics.

```{r}
# HTG
summary <- HTG_calculate_summary_stats(counts_data_tutorial, pattern = "^NC-|^POS-|^GDNA-|^ERCC-")
summary

# RNAseq
summary <- HTG_calculate_summary_stats(counts_data_RNAseq)
summary
```

#### 2.2.3.3 HTG_subset
The `HTG_subset` function is a versatile tool for quickly extracting specific genes or probes from your data, whether it be raw counts or results from differential expression analysis. This functionality is especially useful for focusing on individual genes of interest and facilitating detailed visualization. Additionally, you can request normalization by TPM (Transcripts Per Million) to ensure that the data is appropriately scaled for comparative analysis.

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

* Differential Expression Analysis (DEA)
* Gene Set Enrichment Analysis (GSEA)
* Tumor Microenvironment Analysis (TME)
* Survival Analysis

Since HTG data may have outliers while RNA-seq data typically does not, we will provide two examples for each type of analysis specific to HTG and RNA-seq datasets. 

**NOTE**: The documentation will include various adjustable parameters for each function, allowing you to tailor the analysis to your specific needs.


#### 2.2.4.1 HTG_analysis
All these analysis can be performed using the `HTG_analysis` function.  This function facilitates the execution of differential expression analysis, Gene Set Enrichment Analysis (GSEA), tumor microenvironment analysis, and survival analysis. 

In the two examples provided in this tutorial, the function is configured to perform all analysis. However, you can customize the function to execute only the analysis you require by setting the relevant parameters to `TRUE` or `FALSE`.

**NOTE**: if you don't have a column and is used for an analysis that you will not perform you can put `NULL` on that column

```{r}
# EXAMPLE HTG:
ALL_analysis <- HTG_analysis(
  outliers = outliers_tutorial, 
  pattern = "^NC-|^POS-|^GDNA-|^ERCC-", 
  counts_data = counts_data_tutorial, 
  col_data = AnnotData_tutorial, 
  design_formula = "HPV_status",
  heatmap_columns = c("HPV_status", "Ciclina_D1"), 
  contrast = c("HPV_status", "Positive", "Negative"), 
  pCutoff = 5e-2, 
  variable_01 = "Recurrence_01", 
  time = "Time_to_death_surv", 
  DEA = TRUE,
  remove_outliers = TRUE,
  GSEA = TRUE,
  generate_heatmap = TRUE,
  TME = TRUE,
  survival_analysis = TRUE)


# EXAMPLE RNAseq:
ALL_analysis <- HTG_analysis(
  counts_data = counts_data_RNAseq, 
  col_data = AnnotData_tutorial, 
  design_formula = "HPV_status",
  heatmap_columns = c("HPV_status", "Ciclina_D1"), 
  contrast = c("HPV_status", "Positive", "Negative"), 
  pCutoff = 5e-2, 
  variable_01 = "Recurrence_01", 
  time = "Time_to_death_surv", 
  DEA = TRUE,
  remove_outliers = FALSE,
  GSEA = TRUE,
  generate_heatmap = TRUE,
  TME = TRUE,
  survival_analysis = TRUE)
```

#### 2.2.4.2 HTG_DEA
Each analysis can be performed separately, giving you greater control and flexibility over the process. For example, the `HTG_DEA` function allows you to tailor the differential expression analysis to meet specific needs and datasets. It provides options to apply filters and perform log-fold change (LFC) shrinkage as required, enabling you to customize the analysis according to your objectives.

```{r}
# EXAMPLE FOR HTG: 
results <- HTG_DEA(
  outliers = outliers_tutorial,
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
  counts_data = counts_data_RNAseq,
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
TME_data<- HTG_TME(outliers = outliers_tutorial, pattern =  "^NC-|^POS-|^GDNA-|^ERCC-",
counts_data = counts_data_tutorial , AnnotData= AnnotData_tutorial, design_formula = "HPV_status")

TME_data$EPIC
TME_data$QTI
TME_data$XCELL


# FOR RNAseq
TME_data<- HTG_TME(counts_data = counts_data_RNAseq, AnnotData = AnnotData_tutorial, design_formula = "HPV_status")

TME_data$EPIC
TME_data$QTI
TME_data$XCELL
```

#### 2.2.4.5 HTG_survival
The `HTG_survival` function is designed for performing survival analysis and offers flexibility for various use cases:

* **Custom Genes**: Perform survival analysis on a selected set of genes that you specify.
* **Top 10 Genes from Differential Expression Analysis**: Conduct survival analysis on the top 10 genes identified from your differential expression results.
* **Tumor Microenvironment Analysis Results**: Analyze survival data based on results from tumor microenvironment studies.


```{r}
####   FOR HTG:

# Survial of LCP1 and OMA1 gene
survival_gene_to_use<- HTG_survival(variable_01 = "Recurrence_01",
time = "Time_to_death_surv",
col_data = AnnotData_tutorial,
counts_data = counts_data_tutorial,
res = res_tutorial,
genes_to_use = c("LCP1", "OMA1"),
outliers = outliers_tutorial,
pattern = "^NC-|^POS-|^GDNA-|^ERCC-",
remove_outliers = TRUE)

# Survival from Top10 genes
survival_res<- HTG_survival(variable_01 = "Recurrence_01",
time = "Time_to_death_surv",
col_data = AnnotData_tutorial,
counts_data = counts_data_tutorial,
res = res_tutorial,
genes_to_use = NULL,
outliers = outliers_tutorial,
pattern = "^NC-|^POS-|^GDNA-|^ERCC-",
remove_outliers = TRUE)

# Survival of EPIC TME results.
survival_TME<- HTG_survival(variable_01 = "Recurrence_01",
time = "Time_to_death_surv",
col_data = AnnotData_tutorial,
counts_data = counts_data_tutorial,
res = NULL,
genes_to_use = NULL,
TME = TME_data_tutorial$EPIC,
outliers = outliers_tutorial,
pattern = "^NC-|^POS-|^GDNA-|^ERCC-",
remove_outliers = TRUE)

####   FOR RNAseq: 
                 
# Survial of of CCND1 gene
survival_gene_to_use<- HTG_survival(variable_01 = "Recurrence_01",
time = "Time_to_death_surv",
col_data = AnnotData_tutorial,
counts_data = counts_data_RNAseq,
genes_to_use = "CCND1",
remove_outliers = FALSE)

# Survival from Top10 genes
survival_res<- HTG_survival(variable_01 = "Recurrence_01",
time = "Time_to_death_surv",
col_data = AnnotData_tutorial,
counts_data = counts_data_RNAseq,
res = res_tutorial,
remove_outliers = FALSE)

# Survival of EPIC TME results.
survival_TME<- HTG_survival(variable_01 = "Recurrence_01",
time = "Time_to_death_surv",
col_data = AnnotData_tutorial,
counts_data = counts_data_RNAseq,
res = NULL,
genes_to_use = NULL,
TME = TME_data_tutorial$EPIC,
outliers = outliers_tutorial,
remove_outliers = FALSE)
```
<div style="display: flex; justify-content: space-between;">
    <img src="https://github.com/user-attachments/assets/25dbac67-84eb-4c58-af88-b7e67fdaec33" alt="Image 1" width="300"/>
   <img src="https://github.com/user-attachments/assets/ecf90cb3-9d11-46f7-8b63-cc5c3596902d" alt="ISGlobal Logo" width="300"/>
    <img src="https://github.com/user-attachments/assets/e2680f9a-38e4-4966-bb66-741d2cf58391" alt="Image 2" width="300"/>
</div>
