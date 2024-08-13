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

# 3. QUICK START.

HTGAnalyzer has two functions called HTG_auto and HTG_analysis which provide an easy way to perform all the analyses. Let's start by discussing HTG_auto.

## 3.1 HTG_auto.
This function performs all the analyses mentioned before: Quality Control, Differential Expression Analysis (DEA), Gene Set Enrichment Analysis (GSEA), Tumor Microenvironment Analysis (TME), and Survival Analysis. All the parameters can be modified, but some of them are already set to make it easier. You can refer to the documentation to change them.
```{r}
HTG_auto <- function(counts_file_path,
                     AnnotData_file_path,
                     design_formula,
                     heatmap_columns = NULL,
                     contrast = NULL,
                     variable_01 = NULL,
                     time = NULL)
```
It has more parameters that can be modified, but to keep things simple, let's focus on these:

* `counts_file_path` (Character): Path to the file containing the HTG counts data in Excel format.
* `file_type` (Character): Type of file being imported, either "HTG" or "RNAseq".
* `AnnotData_file_path` (Character): Path to the file containing the annotation data in Excel format.
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
                     "~/AnnotData.xlsx",
                     design_formula = "HPV_status",
                     heatmap_columns = c("HPV_status", "Ciclina_D1"),
                     contrast = c("HPV_status", "Positive", "Negative"),
                     variable_01 = "Recurrence_01",
                     time = "Time_to_death_surv")
```
This will perform a quality control on your data, taking into account the transcriptomic HTG Edge parameters, and will perform a differential expression analysis on HPV_status (Positive vs Negative). It will also generate a heatmap with the columns HPV_status and Ciclina_D1. The results of the differential expression analysis will be used to perform a series of GSEA. Additionally, these results will be used to perform a tumor microenvironment analysis. Finally, it will use the top genes from the differential expression analysis and Recurrence_01 and Time_to_death_surv to perform a survival analysis.


The counts will be imported with the funtion 
```{r}
counts<- HTG_import("path_to_HTG_database.xlsx")
head(counts)
```
After importing, your data will look like this:

|               | 4-B00-03341-B  | 4-B04-23636-B  | 4-B05-03271-B  | 4-B06-00420 | 4-B23-15422 | 4-B23-37783-1 |
|---------------|----------------|----------------|----------------|-------------|-------------|---------------|
| A1BG          | 19             | 24             | 261            | 201         | 14          | 1             |
| A1CF          | 0              | 0              | 429            | 0           | 2           | 0             |
| A2M           | 1472           | 2417           | 4409           | 1990        | 658         | 20            |
| A2ML1         | 0              | 150            | 0              | 46          | 0           | 2             |
| A3GALT2       | 0              | 0              | 0              | 0           | 0           | 0             |

Once you have the counts imported in the desired format, your table will have sample names as column names and genes as row names.
```{r}
colnames(counts)
rownames(counts)
```
You will also find probes at the end:
```{r}
tail(counts)
```

##  2.2 SUBSET
As your data contains probes, you can use the function subset to obtain a data frame with the row names that start with a specific prefix and all the columns (samples).
```{r}
ERCC <- HTG_subset_counts(counts, "ERCC")
NC <- HTG_subset_counts(counts, "NC")
POS <- HTG_subset_counts(counts, "POS")
GDNA <- HTG_subset_counts(counts, "GDNA")
```

## 2.3 FILTERED COUNTS
This function filters counts data to remove rows with specific prefixes: "NC-", "POS-", "GDNA-", and "ERCC-".
```{r}
filtered<- HTG_filterCounts(counts)
```
With the tail function, you can check if all the probes were deleted:
```{r}
tail(counts_filtered)
```

## 2.4 CALCULATE RATIOS
This function calculates ratios based on counts data for different categories such as positive controls and genomic DNA. Results will be stored in a .csv file.
```{r}
ratio<- HTG_calculate_ratios(counts_filtered,POS,NC,GDNA,ERCC)
# check the pathway:
getwd()
head(ratio)
```

## 2.5 QUALITY CONTROL.
This function generates multiple plots to visualize various quality control (QC) metrics. The function includes several thresholds that can be modified based on specific needs, though default values have been tested and are provided. Outliers are marked in red, samples close to the threshold are in yellow, and samples that pass QC are in blue. The quality controls performed are the following:

Density plot of library size.
* QC0 = positive control; threshold = 4 (samples between 3 and 5 are in yellow)
* QC1 = library size; threshold = 7e+06 (samples between 5e+06 and 8e+06 are in yellow)
* QC2 = Negative Control; threshold = 0.045 (samples between 0.035 and 0.05 are in yellow)
* QC3 = Genomic DNA; threshold = 0.02 (samples between 0.015 and 0.025 are in yellow)
* QC4 = ERCC Control; threshold = 5 (samples between 0.015 and 0.03 are in yellow)
* QC5 = Median; threshold = 0.025 (samples between 3 and 7 are in yellow)
* PCA on genes only

```{r}
HTG_plotPCA_probes(counts_filtered)
```

## 2.6 HEATMAP: 
This plot will summarize all the results. In this heatmap, samples that did not pass the control appear in red, and a vector with the samples classified as outliers at least once will be returned.
```{r}
outliers<- HTG_HeatmapQC(ratio,filtered, n_samples = 3)
outliers
```
## 2.7 SUMMARY STATS:
This function calculates summary statistics including minimum, maximum, mean, and median for each column of the input data.
```{r}
HTG_calculate_summary_stats(filtered)
```

## 2.8 PCA Plot with Probes
This performs PCA, identifies outlier samples, and generates plots for PCA results, explained variance, and accumulated variance. You can label the number of samples farthest from the center (if not specified, it will label 3).
```{r}
HTG_plotPCA_genes(counts, 4)
```

## 2.9 PCA ON GENES:
It performs the PCA, identifies outlier samples, and generates plots for PCA results, explained variance, and accumulated variance. You can label the number of samples farthest from the center (if not specified, it will label 3).
```{r}
plotPCA_genes(filtered)
```

# 3. DIFFERENTIAL EXPRESSION ANALYSIS

To perform the Differential Expression analysis you will need a clinical data. 
```{r}
library(readxl)
clinical <- read_excel("path_to_clinical_data.xlsx")
```
Ensure it is a data frame with the Sample ID in rows and also a column with the sample IDs. It is recommended to have no spaces in variable names, which you can change in Excel. It should look like this:
```{r}
head(clinical)
```
|               | id             | Ciclina    | site         | Smoker   |
|---------------|----------------|------------|--------------|----------|
| 4-B00-16661-B | 4-B00-16661-B  | 70         |  lymph node  | yes      |
| 4-B04-28836-B | 4-B04-28836-B  | 70         |  brain       | no       |
| 4-B05-19888-B | 4-B05-19888-B  | 70         | brain        | no       |
| 4-B06-01220   | 4-B06-01220    | 60         | brain        | yes      |
| 4-B16-155432  | 4-B16-155432   | 80         |  lymph node  | no       |
| 4-B16-38476-1 | 4-B16-38476-1  | 40         |  lymph node  | yes      |

As you can see, it can have quantitative and qualitative data.
Three useful functions to better understand the data are:

* str(clinical) = gives the number of rows and columns, the class of each variable, and an example of how it looks.
* summary(clinical) = provides summary statistics for quantitative data.
* table(is.na(clinical$Smoker)) = returns a table with the number of observations that are NA in the specified column (e.g., Smoker).

## 3.1 FROM QUANTITATIVE VARIABLE TO QUALITATIVE VARIABLE
If needed, there is a function to transform quantitative data into qualitative data. If not specified, it will transform the variable into "high" or "low", but this can be changed.
```{r}
table(clinical$Ciclina)
clinical <- HTG_quant_to_qual(clinical, "Ciclina", 34, "high", "low")
head(clinical)
```

## 3.2 DELETE OUTLIERS
Now delete those samples you know are wrong. You can use the function HeatmapQC that also returns the samples.
```{r}
outliers
dim(filtered) 
dim(counts) 
cleaned_counts <- HTG_remove_outliers_COUNTS(filtered, outliers)
annot<-HTG_remove_outliers_ANNOT(clinic,"id",outliers)
dim(cleaned_counts)
dim(annot)
```
## 3.3 SPACE
Variables should not have spaces, so you can replace spaces with underscores using the following code:
```{r}
colnames(annot) <- gsub(" ", "_", colnames(annot))
```
Be careful if there is NA in the variable of interest because it will not work if there is any NA.


## 3.4 CHECK
This step is useful to see if the column of sample IDs has the same names as the columns in our counts data frame.
```{r}
annot <- annot[order(annot$id), ]
cleaned_counts <- cleaned_counts[, order(colnames(cleaned_counts))]
identical(colnames(cleaned_counts), annot$id) #It have to say true.
```

## 3.5 DIFFERENTIAL EXPRESSION ANALYSIS.

This function processes DESeq2 analysis with customizable parameters. It will bring plots with the results to help understand the results as well as the results of the contrast and the distribution of counts on the top ten genes from results.
```{r}
results <- HTG_DESeq(cleaned_counts, annot, "Ciclina",
                               threshold_gene = 200, threshold_subject = 10,
                               heatmap_columns = c("Ciclina", "Smoker"),
                               contrast = c('Ciclina', 'high','low'), pCutoff = 5e-2)
```
# 4. GENE SET ENRICHMENT ANALYSIS

This step will perform and store the results of the gene set enrichment analysis. This is the longest step and will take a while.
```{r}
GSEanalysis<- HTG_GSEAresults(results)
```



Here you will find a tutorial to perform a complete analysis for HTG Edge results with HTGAnalizer.

Having the results of HTG Edge this pakcage will help to perform the Quality Control, the Differential Expression analysis and the Gene Set Enrichment Analysis.
