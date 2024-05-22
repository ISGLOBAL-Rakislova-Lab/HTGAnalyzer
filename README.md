Here you will find a tutorial to perform a complete analysis for HTG Edge results with HTGAnalizer.

Having the results of HTG Edge this pakcage will help to perform the Quality Control, the Differential Expression analysis and the Gene Set Enrichment Analysis.
# INSTALATION:
```{r}
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

library(devtools)
install_github("ISGLOBAL-Rakislova-Lab/HTGAnalizer")
```



# QUALITY CONTROL

## IMPORT
Before starting, make sure your Excel file looks like this. Here we import the data from HTG Edge. This is an Excel file with this structure:

| Sample Name   | 4-B00-03341-B  | 4-B04-23636-B  | 4-B05-03271-B  | 4-B06-00420 | 4-B23-15422 | 4-B23-37783-1 |
|---------------|----------------|----------------|----------------|-------------|-------------|---------------|
| Sample ID     | 1              | 2              | 3              | 4           | 5           | 6             |
| Well          | A1             | B1             | C1             | D1          | E1          | F1            |
| Date Parsed   | 23/04/2022     | 23/04/2022     | 23/04/2022     | 23/04/2022  | 23/04/2022  | 23/04/2022    |
| Total Counts  | 2373           | 6616           | 8720           | 3309        | 1521        | 646           |
| A1BG          | 19             | 24             | 261            | 201         | 14          | 1             |
| A1CF          | 0              | 0              | 429            | 0           | 2           | 0             |
| A2M           | 1472           | 2417           | 4409           | 1990        | 658         | 20            |
| A2ML1         | 0              | 150            | 0              | 46          | 0           | 2             |
| A3GALT2       | 0              | 0              | 0              | 0           | 0           | 0             |
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

##  SUBSET
As your data contains probes, you can use the function subset to obtain a data frame with the row names that start with a specific prefix and all the columns (samples).
```{r}
ERCC <- HTG_subset_counts(counts, "ERCC")
NC <- HTG_subset_counts(counts, "NC")
POS <- HTG_subset_counts(counts, "POS")
GDNA <- HTG_subset_counts(counts, "GDNA")
```

## FILTERED COUNTS
This function filters counts data to remove rows with specific prefixes: "NC-", "POS-", "GDNA-", and "ERCC-".
```{r}
filtered<- HTG_filterCounts(counts)
```
With the tail function, you can check if all the probes were deleted:
```{r}
tail(counts_filtered)
```

## CALCULATE RATIOS
This function calculates ratios based on counts data for different categories such as positive controls and genomic DNA. Results will be stored in a .csv file.
```{r}
ratio<- HTG_calculate_ratios(counts_filtered,POS,NC,GDNA,ERCC)
# check the pathway:
getwd()
head(ratio)
```

## QUALITY CONTROL.
This function generates multiple plots to visualize various quality control (QC) metrics. The function includes several thresholds that can be modified based on specific needs, though default values have been tested and are provided. Outliers are marked in red, samples close to the threshold are in yellow, and samples that pass QC are in blue. The quality controls performed are the following:

Density plot of library size.
* QC0 = positive control; threshold = 4 (samples between 3 and 5 are in yellow)
* QC1 = library size; threshold = 7e+06 (samples between 5e+06 and 8e+06 are in yellow)
* QC2 = Negative Control; threshold = 0.045 (samples between 0.035 and 0.05 are in yellow)
* QC3 = Genomic DNA; threshold = 0.02 (samples between 0.015 and 0.025 are in yellow)
* QC4 = ERCC Control; threshold = 5 (samples between 0.015 and 0.03 are in yellow)
* QC5 = Median; threshold = 0.025 (samples between 3 and 7 are in yellow)

```{r}
HTG_plotPCA_probes(counts_filtered)
```
## SUMMARY STATS:
This function calculates summary statistics including minimum, maximum, mean, and median for each column of the input data.
```{r}
HTG_calculate_summary_stats(a)
```
## HEATMAP: 
This plot will summarize all the results. In this heatmap, samples that did not pass the control appear in red, and a vector with the samples classified as outliers at least once will be returned.
```{r}
outliers<- HTG_HeatmapQC(ratio,filtered, n_samples = 3)
outliers
```

## PCA Plot with Probes
This performs PCA, identifies outlier samples, and generates plots for PCA results, explained variance, and accumulated variance. You can label the number of samples farthest from the center (if not specified, it will label 3).
```{r}
HTG_plotPCA_genes(counts, 4)
```

## PCA ON GENES:
It performs the PCA, identifies outlier samples, and generates plots for PCA results, explained variance, and accumulated variance. You can label the number of samples farthest from the center (if not specified, it will label 3).
```{r}
plotPCA_genes(filtered)
```

# DIFFERENTIAL EXPRESSION ANALYSIS

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

## FROM QUANTITATIVE VARIABLE TO QUALITATIVE VARIABLE
If needed, there is a function to transform quantitative data into qualitative data. If not specified, it will transform the variable into "high" or "low", but this can be changed.
```{r}
table(clinical$Ciclina)
clinical <- HTG_quant_to_qual(clinical, "Ciclina", 34, "high", "low")
head(clinical)
```

## DELETE OUTLIERS
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
## SPACE
Variables should not have spaces, so you can replace spaces with underscores using the following code:
```{r}
colnames(annot) <- gsub(" ", "_", colnames(annot))
```
Be careful if there is NA in the variable of interest because it will not work if there is any NA.


## CHECK
This step is useful to see if the column of sample IDs has the same names as the columns in our counts data frame.
```{r}
annot <- annot[order(annot$id), ]
cleaned_counts <- cleaned_counts[, order(colnames(cleaned_counts))]
identical(colnames(cleaned_counts), annot$id) #It have to say true.
```

## DIFFERENTIAL EXPRESSION ANALYSIS.

This function processes DESeq2 analysis with customizable parameters. It will bring plots with the results to help understand the results as well as the results of the contrast and the distribution of counts on the top ten genes from results.
```{r}
results <- HTG_DESeq(cleaned_counts, annot, "Ciclina",
                               threshold_gene = 200, threshold_subject = 10,
                               heatmap_columns = c("Ciclina", "Smoker"),
                               contrast = c('Ciclina', 'high','low'), pCutoff = 5e-2)
```
# GENE SET ENRICHMENT ANALYSIS

This step will perform and store the results of the gene set enrichment analysis. This is the longest step and will take a while.
```{r}
GSEanalysis<- HTG_GSEAresults(results)
```
