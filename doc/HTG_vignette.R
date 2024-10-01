## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, echo=FALSE--------------------------------------------------------
# Lista de paquetes de Imports
imports <- c(
  "ggplot2", "stats", "apeglm", "data.table", "SummarizedExperiment", "ggrepel", 
  "ggridges", "scales", "ggupset", "pheatmap", "DESeq2", "PoiClaClu", 
  "EnhancedVolcano", "clusterProfiler", "RColorBrewer", "survival", "dplyr", 
  "ggpubr", "tidyr", "maxstat", "readxl", "utils", "graphics", "grDevices", 
  "purrr", "reshape2", "rlang", "tibble"
)

# Lista de paquetes de Suggests
suggests <- c(
  "clusterProfiler", "enrichplot", "ComplexHeatmap", "IOBR", "immunedeconv", 
  "EPIC", "xCell", "org.Hs.eg.db", "knitr", "rmarkdown", "testthat"
)

# Lista de paquetes de Remotes
remotes <- c(
  "omnideconv/immunedeconv", "dviraran/xCell", "GfellerLab/EPIC", 
  "IOBR/IOBR", "kevinblighe/EnhancedVolcano"
)

# Función para cargar un paquete de manera condicional
load_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("El paquete", pkg, "no está instalado."))
  } else {
    library(pkg, character.only = TRUE)
  }
}

# Cargar los paquetes de Imports
for (pkg in imports) {
  load_package(pkg)
}

# Cargar los paquetes de Suggests
for (pkg in suggests) {
  load_package(pkg)
}
library(HTGAnalyzer)


## ----eval=FALSE---------------------------------------------------------------
#  result <- HTG_auto(
#    counts_file_path = "~/HTGAnalyzer/counts.xlsx",  # path to counts data
#    file_type = "HTG",  # Type of file
#    AnnotData_file_path = "~/HTGAnalyzer/annotdata.xlsx",  # path to annotdata.
#    pattern = "^NC-|^POS-|^GDNA-|^ERCC-",  # HTG probes that start with one of these prefixes
#    # this thresholds are stablish for HTG transcriptomic panel.
#    threshold_superior_pos = 5,
#    threshold_line_pos = 4,
#    threshold_inferior_lib = 5e+06,
#    threshold_lib = 7e+06,
#    threshold_superior_nc = 0.05,
#    threshold_line_nc = 0.045,
#    threshold_superior_gdna = 0.025,
#    threshold_line_gdna = 0.02,
#    threshold_superior_ercc = 0.03,
#    threshold_line_ercc = 0.025,
#    threshold_inferior_median = 3,
#    threshold_line_median = 5,
#    save_csv = TRUE,
#    csv_file = "QC_results.csv",
#    design_formula = "HPV_status",
#    percentage_gene = 0.2,
#    threshold_gene = 200,
#    threshold_subject = 10,
#    genes_to_use = c("CCND1", "MMP10", "CTTN"), # Genes used for survival analysis (if you turn NULL it will use the top 10 genes of differential expression. )
#    heatmap_columns = c("HPV_status", "Ciclina_D1"),  # Columns for the heatmap
#    contrast = c("HPV_status", "Positive", "Negative"),  # Contrast for differential expression analysis
#    pCutoff = 5e-2,
#    variable_01 = "Recurrence_01",   # Variable for survival analysis
#    time = "Time_to_death_surv",  # Time for survival analysis
#    # Analysis that can be performed.
#    DEA = TRUE,
#    remove_outliers = TRUE,
#    GSEA = TRUE,
#    generate_heatmap = TRUE,
#    TME = TRUE,
#    survival_analysis = TRUE
#  )
#  

## ----eval=FALSE---------------------------------------------------------------
#  ## LET'S SKIP survival analysis and GSEA
#  
#  result <- HTG_auto(
#    counts_file_path = "~/HTGAnalyzer/counts.xlsx",  # Path to counts data
#    file_type = "HTG",  # Type of file
#    AnnotData_file_path = "~/HTGAnalyzer/annotdata.xlsx",  # Path to annotation data
#    pattern = "^NC-|^POS-|^GDNA-|^ERCC-",  # HTG probes that start with one of these prefixes
#    # These thresholds are specific for the HTG transcriptomic panel
#    threshold_superior_pos = 5,
#    threshold_line_pos = 4,
#    threshold_inferior_lib = 5e+06,
#    threshold_lib = 7e+06,
#    threshold_superior_nc = 0.05,
#    threshold_line_nc = 0.045,
#    threshold_superior_gdna = 0.025,
#    threshold_line_gdna = 0.02,
#    threshold_superior_ercc = 0.03,
#    threshold_line_ercc = 0.025,
#    threshold_inferior_median = 3,
#    threshold_line_median = 5,
#    save_csv = TRUE,
#    csv_file = "QC_results.csv",
#    design_formula = "HPV_status",
#    percentage_gene = 0.2,
#    threshold_gene = 200,
#    threshold_subject = 10,
#    genes_to_use = c("CCND1", "MMP10", "CTTN"),  # Genes used for survival analysis (optional)
#    heatmap_columns = c("HPV_status", "Ciclina_D1"),
#    contrast = c("HPV_status", "Positive", "Negative"),  # Contrast for differential analysis
#    pCutoff = 5e-2,
#    variable_01 = "Recurrence_01",  # Variable for survival analysis
#    time = "Time_to_death_surv",  # Time for survival analysis
#    DEA = TRUE,
#    remove_outliers = TRUE,
#    GSEA = FALSE,
#    generate_heatmap = TRUE,
#    TME = TRUE,
#    survival_analysis = FALSE  # Skip survival analysis by setting this parameter to FALSE
#  )

## ----eval=FALSE---------------------------------------------------------------
#  ## LOAD COUNTS DATA
#  
#  htg_data <- HTG_import_counts("~/HTGAnalyzer/counts.xlsx", "HTG")
#  head(htg_data)
#  
#  ## LOAD ANNOTDATA
#  annotdata<- readxl::read_excel("~/HTGAnalyzer/annotdata.xlsx")
#  head(annotdata)
#  

## -----------------------------------------------------------------------------
head(AnnotData_tutorial)

AnnotData_tutorial <- HTG_quant_to_qual(AnnotData_tutorial,  # Dataset to modify
                                        "Ciclina_D1",  # Column to transform
                                        60,  # Threshold to categorize the data
                                        "high",  # Value for entries greater than 60
                                        "low")  # Value for entries less than 60
#Last column will be the generated one 
head(AnnotData_tutorial)
head(AnnotData_tutorial$Ciclina_D12)


## -----------------------------------------------------------------------------
# Subset rows with the prefix "ERCC"
cts_ERCC <- HTG_subset(counts_data_tutorial, "ERCC")

# Subset rows with the prefix "POS"
cts_POS <- HTG_subset(counts_data_tutorial, "POS")

 # Subset rows with the prefix "NC"
cts_NC <- HTG_subset(counts_data_tutorial, "NC")
 
# Subset rows with the prefix "GDNA"
cts_GDNA <- HTG_subset(counts_data_tutorial, "GDNA")

# Subset rows with the prefix "ZZZ3"
cts_ZZZ3 <- HTG_subset(counts_data_tutorial, "ZZZ3")

# Subset rows with the prefix "ZZZ3" and normalize the data
cts_ZZZ3 <- HTG_subset(counts_data_tutorial, "ZZZ3", normalize = TRUE)

resAAAS<- HTG_subset(res_tutorial, "AAAS")

## -----------------------------------------------------------------------------
outliers <- HTG_QC(counts_data_tutorial, pattern = "^NC-|^POS-|^GDNA-|^ERCC-", 
                             threshold_superior_pos = 5,
                             threshold_line_pos = 4,
                             threshold_inferior_lib = 5e+06,
                             threshold_lib = 7e+06,
                             threshold_superior_nc = 0.05,
                             threshold_line_nc = 0.045,
                             threshold_superior_gdna = 0.025,
                             threshold_line_gdna = 0.02,
                             threshold_superior_ercc = 0.03,
                             threshold_line_ercc = 0.025,
                             threshold_inferior_median = 3,
                             threshold_line_median = 5,
                             save_csv = TRUE,
                             csv_file = "QC_results.csv")


## -----------------------------------------------------------------------------
summary <- HTG_calculate_summary_stats(counts_data_tutorial, pattern = "^NC-|^POS-|^GDNA-|^ERCC-")

## -----------------------------------------------------------------------------
ALL_done <- HTG_analysis(outliers = outliers_tutorial, 
                         pattern = "^NC-|^POS-|^GDNA-|^ERCC-", 
                         counts_data = counts_data_tutorial, 
                         col_data = AnnotData_tutorial, 
                         design_formula = "HPV_status", 
                         percentage_gene = 0.2, 
                         threshold_gene = 200, 
                         threshold_subject = 10,
                         genes_to_use = c("CCND1", "MMP10", "CTTN"),
                         heatmap_columns = c("HPV_status","Ciclina_D1"),
                         contrast = c("HPV_status", "Positive", "Negative"), 
                         pCutoff = 5e-2, 
                         variable_01 = "Recurrence_01",
                         time = "Time_to_death_surv", 
                         DEA = TRUE, 
                         remove_outliers = TRUE,
                         GSEA = FALSE, 
                         generate_heatmap = TRUE, 
                         TME = TRUE, 
                         survival_analysis = FALSE)

## -----------------------------------------------------------------------------
### DEA:
results <- HTG_DEA(outliers= outliers_tutorial, 
                   counts_data= counts_data_tutorial, 
                   AnnotData_tutorial, 
                   design_formula = "HPV_status", 
                   heatmap_columns = c("HPV_status", "Recurrence"),
                   contrast = c("HPV_status", "Positive", "Negative"), 
                   pattern = "^NC-|^POS-|^GDNA-|^ERCC-", 
                   remove_outliers = TRUE,
                   percentage_gene = 0.2, 
                   threshold_gene = 200,
                   threshold_subject = 10,
                   pCutoff = 5e-2, 
                   apply_filtering = TRUE, 
                   apply_lfc_shrinkage = TRUE, 
                   extract_shrinkage = FALSE)

head(results)

## -----------------------------------------------------------------------------
HTG_GSEA(res_tutorial)

## -----------------------------------------------------------------------------

TME_results<- HTG_TME(outliers= outliers_tutorial, 
                      pattern= "^NC-|^POS-|^GDNA-|^ERCC-",
                      counts_data= counts_data_tutorial, 
                      AnnotData= AnnotData_tutorial, 
                      design_formula= "HPV_status", 
                      remove_outliers = TRUE, 
                      DEA = NULL)


head(TME_results$EPIC)
head(TME_results$QTI)
head(TME_results$XCELL)

## -----------------------------------------------------------------------------
############# survival using res
 HTG_survival(variable_01 = "Recurrence_01",
              time = "Time_to_death_surv",
              col_data = AnnotData_tutorial,
              counts_data = counts_data_tutorial,
              res = res_tutorial, #DEA results. 
              genes_to_use = NULL,
              outliers = outliers_tutorial,
              pattern = "^NC-|^POS-|^GDNA-|^ERCC-",
              remove_outliers = TRUE)

###### survival using genes to use
 HTG_survival(variable_01 = "Recurrence_01",
              time = "Time_to_death_surv",
              col_data = AnnotData_tutorial,
              counts_data = counts_data_tutorial,
              res = NULL,
              genes_to_use = c("LCP1", "OMA1"),
              outliers = outliers_tutorial,
              pattern = "^NC-|^POS-|^GDNA-|^ERCC-",
              remove_outliers = TRUE)
 

###### survival using TME data.
HTG_survival(variable_01 = "Recurrence_01",
              time = "Time_to_death_surv",
              col_data = AnnotData_tutorial,
              counts_data = counts_data_tutorial,
              res = NULL,
              genes_to_use = NULL,
              TME = TME_data_tutorial$EPIC,
              outliers = outliers,
              pattern = "^NC-|^POS-|^GDNA-|^ERCC-",
              remove_outliers = TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  result <- HTG_auto(
#    counts_file_path = "~/HTGAnalyzer/counts_RNAseq.xlsx",  # path to counts data
#    file_type = "RNAseq",  # Tipe of file
#    AnnotData_file_path = "~/HTGAnalyzer/annotdata.xlsx",  # path to annotdata.
#    pattern = NULL,
#    save_csv = TRUE,
#    csv_file = "QC_results.csv",
#    design_formula = "HPV_status",
#    percentage_gene = 0.2,
#    threshold_gene = 200,
#    threshold_subject = 10,
#    genes_to_use = c("CCND1", "MMP10", "CTTN"), # Genes used for survival analysis (optional)
#    heatmap_columns = c("HPV_status", "Ciclina_D1"),  # Columns for the heatmap
#    contrast = c("HPV_status", "Positive", "Negative"),  # Contrast for differential analysis
#    pCutoff = 5e-2,
#    variable_01 = "Recurrence_01",   # Variable for analysis
#    time = "Time_to_death_surv",  # Time for survival analysis
#    DEA = TRUE,
#    remove_outliers = TRUE,
#    GSEA = TRUE,
#    generate_heatmap = TRUE,
#    TME = TRUE,
#    survival_analysis = TRUE
#  )
#  

## ----eval=FALSE---------------------------------------------------------------
#  result <- HTG_auto(
#    counts_file_path = "~/HTGAnalyzer/counts_RNAseq.xlsx",  # path to counts data
#    file_type = "RNAseq",  # Tipe of file
#    AnnotData_file_path = "~/HTGAnalyzer/annotdata_RNAseq.xlsx",  # path to annotdata.
#    pattern = NULL,
#    save_csv = TRUE,
#    csv_file = "QC_results.csv",
#    design_formula = "HPV_status",
#    percentage_gene = 0.2,
#    threshold_gene = 200,
#    threshold_subject = 10,
#    genes_to_use = c("CCND1", "MMP10", "CTTN"),  # Genes used for survival analysis (optional)
#    heatmap_columns = c("HPV_status", "Ciclina_D1"),
#    contrast = c("HPV_status", "Positive", "Negative"),  # Contrast for differential analysis
#    pCutoff = 5e-2,
#    variable_01 = "Recurrence_01",  # Variable for survival analysis
#    time = "Time_to_death_surv",  # Time for survival analysis
#    DEA = TRUE,
#    remove_outliers = TRUE,
#    GSEA = FALSE,
#    generate_heatmap = TRUE,
#    TME = TRUE,
#    survival_analysis = FALSE  # Skip survival analysis by setting this parameter to FALSE
#  )

## ----eval=FALSE---------------------------------------------------------------
#  ## LOAD COUNTS DATA
#  
#  htg_data <- HTG_import_counts("~/HTGAnalyzer/counts.xlsx", "RNAseq")
#  head(htg_data)
#  
#  ## LOAD ANNOTDATA
#  annotdata<- readxl::read_excel("~/HTGAnalyzer/annotdata.xlsx")
#  head(annotdata)
#  

## -----------------------------------------------------------------------------
head(AnnotData_tutorial)

AnnotData_tutorial <- HTG_quant_to_qual(AnnotData_tutorial,  # Dataset to modify
                                        "Ciclina_D1",  # Column to transform
                                        60,  # Threshold to categorize the data
                                        "high",  # Value for entries greater than 60
                                        "low")  # Value for entries less than 60
head(AnnotData_tutorial)


## -----------------------------------------------------------------------------
 # Subset rows with the prefix "ERCC"
 cts_ERCC <- HTG_subset(counts_data_tutorial, "ERCC")

 # Subset rows with the prefix "POS"
 cts_POS <- HTG_subset(counts_data_tutorial, "POS")

 # Subset rows with the prefix "NC"
 cts_NC <- HTG_subset(counts_data_tutorial, "NC")
 
# Subset rows with the prefix "GDNA"
 cts_GDNA <- HTG_subset(counts_data_tutorial, "GDNA")

 # Subset rows with the prefix "ZZZ3"
 cts_ZZZ3 <- HTG_subset(counts_data_tutorial, "ZZZ3")

  # Subset rows with the prefix "ZZZ3" and normalize the data
 cts_ZZZ3 <- HTG_subset(counts_data_tutorial, "ZZZ3", normalize = TRUE)

 resAAAS<- HTG_subset(res_tutorial, "AAAS")

## -----------------------------------------------------------------------------
summary <- HTG_calculate_summary_stats(counts_data_tutorial, pattern = NULL)

## -----------------------------------------------------------------------------
ALL_done <- HTG_analysis(outliers = NULL, # RNAseq doen't have outliers. 
                         pattern = NULL, 
                         counts_data = counts_data_tutorial, 
                         col_data = AnnotData_tutorial, 
                         design_formula = "HPV_status", 
                         percentage_gene = 0.2, 
                         threshold_gene = 200, 
                         threshold_subject = 10,
                         genes_to_use = c("CCND1", "MMP10", "CTTN"),
                         heatmap_columns = c("HPV_status","Ciclina_D1"),
                         contrast = c("HPV_status", "Positive", "Negative"), 
                         pCutoff = 5e-2, 
                         variable_01 = "Recurrence_01",
                         time = "Time_to_death_surv", 
                         DEA = TRUE, 
                         remove_outliers = FALSE, #set false because don't have outliers. 
                         GSEA = FALSE, 
                         generate_heatmap = TRUE, 
                         TME = TRUE, 
                         survival_analysis = FALSE)

## -----------------------------------------------------------------------------
### DEA:
results <- HTG_DEA(outliers= NULL, 
                   counts_data= counts_data_tutorial, 
                   AnnotData_tutorial, 
                   design_formula = "HPV_status", 
                   heatmap_columns = c("HPV_status", "Recurrence"),
                   contrast = c("HPV_status", "Positive", "Negative"), 
                   pattern = NULL, 
                   remove_outliers = TRUE,
                   percentage_gene = 0.2, 
                   threshold_gene = 200,
                   threshold_subject = 10,
                   pCutoff = 5e-2, 
                   apply_filtering = TRUE, 
                   apply_lfc_shrinkage = TRUE, 
                   extract_shrinkage = FALSE)

head(results)

## -----------------------------------------------------------------------------
HTG_GSEA(res_tutorial)

## -----------------------------------------------------------------------------

TME_results<- HTG_TME(outliers= NULL, 
                      pattern= NULL,
                      counts_data= counts_data_tutorial, 
                      AnnotData= AnnotData_tutorial, 
                      design_formula= "HPV_status", 
                      remove_outliers = FALSE, 
                      DEA = NULL)


head(TME_results$EPIC)
head(TME_results$QTI)
head(TME_results$XCELL)

## -----------------------------------------------------------------------------
###### survival using res
 HTG_survival(variable_01 = "Recurrence_01",
              time = "Time_to_death_surv",
              col_data = AnnotData_tutorial,
              counts_data = counts_data_tutorial,
              res = res_tutorial, #DEA results. 
              genes_to_use = NULL,
              outliers = NULL,
              pattern = NULL,
              remove_outliers = FALSE)

###### survival using genes to use
 HTG_survival(variable_01 = "Recurrence_01",
              time = "Time_to_death_surv",
              col_data = AnnotData_tutorial,
              counts_data = counts_data_tutorial,
              res = NULL,
              genes_to_use = c("LCP1", "OMA1"),
              outliers = NULL,
              pattern = NULL,
              remove_outliers = FALSE)
 

###### survival using TME data.
HTG_survival(variable_01 = "Recurrence_01",
              time = "Time_to_death_surv",
              col_data = AnnotData_tutorial,
              counts_data = counts_data_tutorial,
              res = NULL,
              genes_to_use = NULL,
              TME = TME_data_tutorial$EPIC,
              outliers = NULL,
              pattern = NULL,
              remove_outliers = FALSE)

