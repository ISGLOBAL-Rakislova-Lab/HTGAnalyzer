Here you will find a tutorial to perform a complete analysis for HTG Edge results and RNAseq with HTGAnalyzer.

Having the results of HTG Edge Reveal or RNAseq, this package will help you perform Quality Control, Differential Expression Analysis (DEA), Gene Set Enrichment Analysis (GSEA), Tumor Microenvironment Analysis (TM), and Survival Analysis.

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
Once installed, you can start with the tutorial. Remember that this package was initially designed for people who have HTG Edge to analyze. However, with the closure of the company, previous quality control features were revealed, but not for the transcriptomic panel, whose quality control was not well established. For this reason, this package is prepared to perform an easy analysis for people who are not into bioinformatics. It has an automatic function that will perform all the analyses and other functions for specific analyses, giving you more freedom to modify.

We have a main function, HTG_auto, which performs all the mentioned analyses. Additionally, we have individual functions such as HTG_DEA, HTG_survival, HTG_QC, and HTG_TME, which provide more flexibility for modification.

## 2.1. QUICK START.

HTGAnalyzer has a functions called HTG_auto which provide an easy way to perform all the analyses. Let's start by discussing HTG_auto.

### 2.1.1 HTG_auto.
The `HTG_auto` function in HTGAnalyzer automates a comprehensive analysis pipeline for HTG data. It integrates various analyses including Quality Control (QC), Differential Expression Analysis (DEA), Gene Set Enrichment Analysis (GSEA), Tumor Microenvironment Analysis (TME), and Survival Analysis. This function is designed to simplify the process, with default settings for ease of use, while also providing flexibility to modify parameters according to specific needs.
```{r}
HTG_auto <- function(counts_file_path,
                     file_type,  
                     AnnotData_file_path,
                     design_formula,
                     QC = TRUE,
                     heatmap_columns = NULL,
                     contrast = NULL,
                     variable_01 = NULL,
                     time = NULL)
```
It has more parameters that can be modified, but to keep things simple, let's focus on these:

* `counts_file_path` (Character): Path to the file containing the HTG counts data in Excel format.
* `file_type` (Character): Type of file being imported, either "HTG" or "RNAseq".
* `AnnotData_file_path` (Character): Path to the file containing the annotation data in Excel format.
* `QC` Indicates whether to perform quality control on the data. Default is TRUE. If set to FALSE, the QC step is skipped (for RNAseq data).
* `design_formula` (Character): The design formula for DESeq2 analysis, specified as a string without the tilde (~).
* `heatmap_columns` (Character vector): Specifies columns to be used for annotations in the heatmap.
* `contrast` (Character vector): Specifies the contrast for differential expression analysis. Default is NULL.
* `variable_01` (Character): Name of the survival event variable (e.g., "Recurrence_01"). Required for survival analysis.
* `time` (Character): Name of the time variable (e.g., "Time_to_death_surv"). Required for survival analysis.
  
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
HTG_auto <- function("~/HPV_counts.xlsx",
                     file_type = "HTG",
                     "~/AnnotData.xlsx",
                     design_formula = "HPV_status",
                     QC = TRUE,
                     heatmap_columns = c("HPV_status", "Ciclina_D1"),
                     contrast = c("HPV_status", "Positive", "Negative"),
                     variable_01 = "Recurrence_01",
                     time = "Time_to_death_surv")
```
Certainly! Here’s a narrative description based on the provided text:

---

In the case of HTG data, setting `QC = TRUE` is crucial as it ensures that the quality control process is carried out, allowing the detection and handling of outliers. This step is essential for validating the integrity of the data before any further analysis. With quality control enabled, the function will proceed to perform a Differential Expression Analysis (DEA) based on the `HPV_status`, comparing positive and negative samples. This comparison will reveal which genes are significantly differentially expressed between the two HPV status groups.

Following the DEA, a heatmap will be generated using the columns `HPV_status` and `Ciclina_D1`. This visualization helps in understanding how well the data segregates based on these variables. By examining the heatmap, you can gain insights into which variables are most effective in distinguishing between the groups, providing a clearer picture of the data’s structure.

The results from the DEA will then be utilized to perform Gene Set Enrichment Analysis (GSEA). GSEA will identify whether specific gene sets are significantly enriched in the differentially expressed genes, offering deeper biological insights into the underlying mechanisms driving the differences observed between HPV-positive and HPV-negative samples. This can help in understanding the biological pathways or processes that are activated or repressed in response to HPV status.

Additionally, the findings from the DEA will be applied to a Tumor Microenvironment (TME) analysis. This analysis can be particularly valuable as it provides insights into the interaction between tumor cells and their surrounding environment. Understanding these interactions can shed light on how the tumor microenvironment influences gene expression and contributes to the overall tumor biology, potentially revealing new targets for therapeutic intervention.

Finally, the top genes identified from the differential expression analysis, along with the survival-related variables `Recurrence_01` and `Time_to_death_surv`, will be used to conduct a survival analysis. This step is critical as it evaluates the impact of gene expression on patient outcomes, providing valuable prognostic information. By correlating gene expression profiles with survival data, you can identify genes associated with better or worse prognosis, which can be crucial for developing personalized treatment strategies and improving patient management.

Keep in mind that this example is tailored for HTG data, where the quality control step is specifically programmed to accommodate the characteristics of transcriptomic panels. However, if you are working with RNAseq data, you will need to adjust some parameters accordingly. The quality control thresholds and other settings should be customized based on the specifics of the RNAseq data to ensure accurate and reliable analysis results.


## 2.2. PRINCIPAL FUNCTIONS AND SECUNDARY FUNCTIONS 
Just in case you need more control of your data or you want to perform other analysis. HTGAnalyzer have principal and secundary fucntions to help reach your demands.

PRINCIPAL: 
* `HTG_import_counts`:Import counts data from an Excel file. It can be either HTG excel or RNAseq. The first row in the Excel file must contain the column headers, with "id" as the first column header followed by the names of each sample. Please note that you might need to modify the Excel file to ensure it is in the correct format for importing data into R. Ensure column names are free of special characters (e.g., spaces, (,), ?, `, ^, ., -, *, and others) to avoid import issues. The program will change the spaces into `_`
* `HTG_QC`: This function performs comprehensive quality control (QC) for HTG EdgeSeq data, including checks on positive values, library size, and various thresholds (e.g., negative control, genomic DNA, ERCC). It generates a data frame with summary statistics and sample ratios, highlights outliers with an optional heatmap, and saves plots in the working directory.
* `HTG_analysis`: This function conducts a comprehensive analysis pipeline including DESeq2 differential expression analysis (DEA), Gene Set Enrichment Analysis (GSEA), TME analysis, and survival analysis. The pipeline supports optional steps for generating volcano plots and heatmaps. The function is suitable for both HTG and RNA-seq data.

SECUNDARY:
* `HTG_subset`: This function subsets a data frame based on a specified prefix in the row names. It allows you to extract specific rows that match a given pattern and optionally normalizes the data using TPM (Transcripts Per Million) before subsetting. The function also provides the dimensions of the resulting data frame.
* `HTG_plotPCA`: This function generates a PCA plot specifically for genes, along with plots showing explained variability and accumulated explained variability. It highlights the samples that are the most distant from the center of the PCA plot.
* `quant_to_qual`: This function converts quantitative columns in a data frame into qualitative ones based on a specified threshold.

### 2.2.1 IN-DEPTH GUIDE

For this part of the tutorial, we will combine principal and secondary functions to explore their capabilities. 

#### 2.2.1.1 DATA IMPORT.

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
head(AnnotData
```
We recomend you to avoir special characters on excel files (e.g., spaces, (,), ?, `, ^, ., -, *, and others)  to avoid analysis issues.

#### 2.2.1.2 QUALITY CONTROL.
This function performs various quality control (QC) checks tailored for the HTG EdgeSeq transcriptomic panel, though thresholds can be adjusted as needed. QC checks include:

* QC0: % of positive values > 4%.
* QC1: Library size > 7e+06.
* QC2: Negative control threshold > 0.045.
* QC3: Genomic DNA threshold > 0.02.
* QC4: ERCC threshold > 0.025.

The function generates a data frame (which can be saved as a .csv file, with a preview shown in the R console) with:
* Sum of each probe for each sample (total genes, positive, negative, gdna, ercc).
* Ratios for each sample.
* Sample sizes.
* A PCA column indicating samples furthest from the center.
Additionally, a statistical .csv is generated with columns for Min, Max, Mean, Median, Mode, SD, Variance, Range, Q1, Q3, IQR, Skewness, Kurtosis, Missing, and CV.

The function includes a heatmap to highlight potential outliers, which will be saved in a vector. Additionally, it performs a PCA analysis that identifies and highlights the specified number of samples that are furthest from the center. The PCA results include plots showing the explained variability and the cumulative explained variability.

```{r}
outliers<- HTG_QC(counts_data)
outliers
```
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
# PCA FOR RNAseq
PCA_RNAseq<- HTG_plotPCA(counts_data, n_samples = 0,pattern = NULL)
PCA_RNAseq
summary <- HTG_calculate_summary_stats(counts_data)
summary

# PCA FOR HTG
PCA_HTG<- HTG_plotPCA(counts_data, n_samples = 0,pattern = "^NC-|^POS-|^GDNA-|^ERCC-")
PCA_HTG
summary <- HTG_calculate_summary_stats(counts_data, pattern = "^NC-|^POS-|^GDNA-|^ERCC-")
summary
```

The `HTG_subset` function is a versatile tool for quickly extracting specific genes or probes from your data, whether it be raw counts or results from differential expression analysis. This functionality is especially useful for focusing on individual genes of interest and facilitating detailed visualization.

EXAMPLES:

```r
# Subset ERCC probes from counts data
  ERCC <- HTG_subset(counts_data, "ERCC")

# Extracting a Specific Gene (AAAS) from Counts Data
# Subset AAAS gene from counts data without normalization
counts_AAAS <- HTG_subset(counts_data, "AAAS")

# Subset AAAS gene from counts data with normalization
counts_AAAS_normalize <- HTG_subset(counts_data, "AAAS", normalize = TRUE)

# Extracting a Specific Gene (AAAS) from Differential Expression Results
# Subset AAAS gene from differential expression results
res_AAAS <- HTG_subset(res, "AAAS")
```

#### 2.2.1.3 ANALYSIS
Once the QC is performed you can also 4 diferent analysis:
* differential expression analysis
* GSEA
* Tumor microenviroment
* survival analysis.

All this analysis can be perfomed with HTG_analysis
```{r}

a <- HTG_analysis(
  outliers = outliers, 
  pattern = "^NC-|^POS-|^GDNA-|^ERCC-", 
  counts_data = counts_data, 
  col_data = clinical, 
  design_formula = "HPV_status",
  percentage_gene = 0.2, 
  percentage_zero = 0.2, 
  threshold_gene = 200, 
  threshold_subject = 10, 
  top_genes = c("CCND1", "MMP10", "CTTN"), 
  heatmap_columns = c("HPV_status", "Ciclina_D1"), 
  contrast = c("HPV_status", "Positive", "Negative"), 
  pCutoff = 5e-2, 
  variable_01 = "Recurrence_01", 
  time = "Time_to_death_surv", 
  dds = TRUE, 
  generate_volcano = FALSE, 
  remove_outliers = TRUE, 
  GSEA = TRUE, 
  generate_heatmap = TRUE, 
  TME = FALSE, 
  survival_analysis = FALSE, 
  grupos = NULL
)

```

