#' HTG_QC: Quality Control for HTG EdgeSeq Data
#'
#' @description
#' This function performs various quality control (QC) checks for HTG EdgeSeq transcriptomic panel data. The QC checks include:
#'
#' QC0: Percentage of positive values < than 4%; QC1: Library size > than 7e+06; QC2: Negative control threshold < than 0.045; QC3: Genomic DNA threshold < than 0.02; QC4: ERCC threshold < than 0.025; Median:threshold > 5.
#'
#' In addition to these plots, highlighting potential outlier samples. The function also creates a data frame including the sum of each probe for each sample (total genes, positive, negative, gdna, and ercc)
#' the ratio for each sample and the size of each sample. Additionally, a statistical .csv is generated with columns for Min, Max, Mean, Median, Mode, SD, Variance, Range, Q1, Q3, IQR, Skewness, Kurtosis, Missing, and CV.
#'
#' This function also includes an optional heatmap to highlight potential outlier samples, which will be saved in a vector.
#'
#' The plots will saved in the current working directory.
#'
#' @param counts_data A data frame containing the HTG count data. The data must include probes that start with "^NC-|^POS-|^GDNA-|^ERCC-" for the function to work correctly.
#' @param pattern A regular expression pattern to identify control probes in the count data. For HTG data, this could be "^NC-|^POS-|^GDNA-|^ERCC-". If NULL, the pattern will not be applied.
#' @param threshold_superior_pos Threshold for upper limit of positive control ratio.
#' @param threshold_line_pos Threshold line for positive control ratio.
#' @param threshold_inferior_lib Threshold for lower limit of library size.
#' @param threshold_lib Threshold line for library size.
#' @param threshold_superior_nc Threshold for upper limit of negative control ratio.
#' @param threshold_line_nc Threshold line for negative control ratio.
#' @param threshold_superior_gdna Threshold for upper limit of genomic DNA ratio.
#' @param threshold_line_gdna Threshold line for genomic DNA ratio.
#' @param threshold_superior_ercc Threshold for upper limit of ERCC control ratio.
#' @param threshold_line_ercc Threshold line for ERCC control ratio.
#' @param threshold_inferior_median Threshold for lower limit of median ratio.
#' @param threshold_line_median Threshold line for median ratio.
#' @param save_csv Logical, whether to save the ratios as a CSV file. Default is FALSE.
#' @param csv_file The name of the CSV file to save the ratios if save_csv is TRUE. Default is "QC_results.csv".
#'
#' @return This function generates multiple plots displaying various QC metrics, including a violin plots, and saves an Excel file with all the ratios. Additionally, it identifies and returns the most probable outliers based on the QC analysis.
#'
#' @export
#'
#'
#' @examples
#' # Run the function with example data
#' HTG_QC(counts_data = counts_data_tutorial, pattern = "^NC-|^POS-|^GDNA-|^ERCC-", save_csv = TRUE)
#' qc_results <- HTG_QC(
#'  counts_data = counts_data_tutorial,
#'  pattern = "^NC-|^POS-|^GDNA-|^ERCC-",   # control probe patterns
#'  threshold_superior_pos = 5,             # maximum allowed value for POS controls
#'  threshold_line_pos = 4,                 # warning limit for POS controls
#'  threshold_inferior_lib = 5e+06,         # minimum acceptable library size
#'  threshold_lib = 7e+06,                  # warning limit for library size
#'  threshold_superior_nc = 0.05,           # maximum allowed for NC (negative controls)
#'  threshold_line_nc = 0.045,              # warning limit for NC
#'  threshold_superior_gdna = 0.025,        # maximum allowed for gDNA
#'  threshold_line_gdna = 0.02,             # warning limit for gDNA
#'  threshold_superior_ercc = 0.03,         # maximum allowed for ERCC
#'  threshold_line_ercc = 0.025,            # warning limit for ERCC
#'  threshold_inferior_median = 3,          # minimum acceptable for counts median
#'  threshold_line_median = 5,              # warning limit for counts median
#'  save_csv = TRUE,                        # save QC results to CSV
#'  csv_file = "QC_results_tutorial.csv"    # output file name
#' )
#'
#' @name HTG_QC
#'
utils::globalVariables(c("PC1", "PC2", "Tag", "label", "Componente", "Porcentaje", "pc", "Sample", "LogTPM", "variable", "value"))
HTG_QC <- function(counts_data, pattern = "^NC-|^POS-|^GDNA-|^ERCC-",
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
                             csv_file = "QC_results.csv") {

  # Filter counts_data data
  cat("\033[33mINITIATING DATA FILTERING...\033[0m\n")
  counts_filtered <- subset(counts_data, !grepl(pattern, rownames(counts_data)))
  min_values <- apply(counts_filtered, 2, min)
  max_values <- apply(counts_filtered, 2, max)
  mean_values <- apply(counts_filtered, 2, mean)
  median_values <- apply(counts_filtered, 2, median)
  mode_values <- apply(counts_filtered, 2, function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  })
  sd_values <- apply(counts_filtered, 2, sd)
  var_values <- apply(counts_filtered, 2, var)
  range_values <- apply(counts_filtered, 2, function(x) max(x) - min(x))
  quartile_1 <- apply(counts_filtered, 2, function(x) quantile(x, 0.25))
  quartile_3 <- apply(counts_filtered, 2, function(x) quantile(x, 0.75))
  iqr_values <- quartile_3 - quartile_1
  skewness_values <- apply(counts_filtered, 2, function(x) {
    n <- length(x)
    mean_x <- mean(x)
    sd_x <- sd(x)
    sum((x - mean_x)^3) / ((n - 1) * (sd_x^3))
  })
  kurtosis_values <- apply(counts_filtered, 2, function(x) {
    n <- length(x)
    mean_x <- mean(x)
    sd_x <- sd(x)
    sum((x - mean_x)^4) / ((n - 1) * (sd_x^4)) - 3
  })
  missing_values <- apply(counts_filtered, 2, function(x) sum(is.na(x)))
  cv_values <- sd_values / mean_values

  summary_stats <- data.frame(
    Min = min_values,
    Max = max_values,
    Mean = mean_values,
    Median = median_values,
    Mode = mode_values,
    SD = sd_values,
    Variance = var_values,
    Range = range_values,
    Q1 = quartile_1,
    Q3 = quartile_3,
    IQR = iqr_values,
    Skewness = skewness_values,
    Kurtosis = kurtosis_values,
    Missing = missing_values,
    CV = cv_values
  )
  summary_stats$ID <- rownames(summary_stats)
  write.csv(summary_stats, file = "summary_stats.csv", row.names = FALSE)
  cat("\033[32mSummary statistics saved as 'summary_stats.csv'\033[0m\n")

  cat("\033[33mINITIATING QC PLOTS...\033[0m\n")
  # Subsets
  cts_ERCC <- as.data.frame(subset(counts_data, grepl("^ERCC-", rownames(counts_data))))
  cts_NC <- as.data.frame(subset(counts_data, grepl("^NC-", rownames(counts_data))))
  cts_POS <- as.data.frame(subset(counts_data, grepl("^POS-", rownames(counts_data))))
  cts_GDNA <- as.data.frame(subset(counts_data, grepl("^GDNA-", rownames(counts_data))))

  # Ratios
  total_gens <- colSums(counts_filtered)
  total_POS <- colSums(cts_POS)
  total_NC <- colSums(cts_NC)
  total_GDNA <- colSums(cts_GDNA)
  total_ERCC <- colSums(cts_ERCC)

  # Calculate ratios
  ratios <- data.frame(
    total_POS = total_POS,
    total_GDNA = total_GDNA,
    total_gens = total_gens,
    total_NC = total_NC,
    total_ERCC = total_ERCC
  )
  ratios$`pos/gens` <- (ratios$total_POS / ratios$total_gens) * 100
  ratios$`gdna/gens` <- (ratios$total_GDNA / ratios$total_gens) * 100
  ratios$`nc/gens` <- (ratios$total_NC / ratios$total_gens) * 100
  ratios$`ERCC/gens` <- (ratios$total_ERCC / ratios$total_gens) * 100

  # Add median column
  summary_stats <- HTG_calculate_summary_stats(counts_filtered, pattern= pattern)
  ratios$median <- summary_stats$Median
  ratiosb<-ratios
  ratiosb$min<- summary_stats$Min
  ratiosb$max<- summary_stats$Max
  ratiosb$mean<- summary_stats$Mean

  # Add sample names as a factor column
  ratios$samples <- factor(rownames(ratios))

  # Calcula los valores de las columnas Q0 a Q5 usando los umbrales definidos
  ratiosb$Q0 <- ifelse(ratiosb$`pos/gens` < threshold_line_pos, "PASS", "FAIL")
  ratiosb$Q1 <- ifelse(ratiosb$total_gens > threshold_lib, "PASS", "FAIL")
  ratiosb$Q2 <- ifelse(ratiosb$total_NC > threshold_line_nc, "PASS", "FAIL")
  ratiosb$Q3 <- ifelse(ratiosb$`gdna/gens` < threshold_line_gdna, "PASS", "FAIL")
  ratiosb$Q4 <- ifelse(ratiosb$`ERCC/gens` < threshold_line_ercc, "PASS", "FAIL")
  ratiosb$Q5 <- ifelse(ratiosb$median > threshold_line_median, "PASS", "FAIL")

  # Agrega los rownames como una columna 'ID'
  ratiosb$ID <- rownames(ratiosb)

  # Extrae el número de paciente y ordena por este número
  ratiosb$ID_num <- as.numeric(gsub("patient_", "", ratiosb$ID))  # Extrae el número de ID
  ratiosb <- ratiosb[order(ratiosb$ID_num), ]  # Ordena por la columna ID_num
  ratiosb$ID_num <- NULL  # Elimina la columna auxiliar

  # Elimina los rownames originales y reorganiza las columnas según el orden deseado
  rownames(ratiosb) <- NULL
  ratiosb <- ratiosb[, c("ID", "total_POS", "total_NC", "total_GDNA", "total_ERCC",
                       "pos/gens", "Q0", "total_gens", "Q1",
                       "nc/gens", "Q2", "gdna/gens", "Q3",
                       "ERCC/gens", "Q4", "median", "Q5",
                       "min", "max", "mean")]
  ratiosb$QC_status <- ifelse(
    apply(ratiosb[, c("Q0", "Q1", "Q2", "Q3", "Q4", "Q5")], 1, function(x) any(x == "FAIL")),
    "FAIL",
    "PASS"
  )

  # Verifica el resultado final
  print(ratiosb)

  # Optionally save as CSV
  if (save_csv) {
    write.csv(ratiosb, csv_file, row.names = FALSE)
    cat("\033[32mQC DATA SAVED AS '", csv_file, "'\033[0m\n")
  }
###
  library_size <- colSums(counts_filtered)

  # Create dataframe of library size
  lib_s2 <- data.frame(Sample = colnames(counts_filtered), Size = library_size)
  ratios_heat <- as.matrix(ratios)

  # Add a fourth column to ratios_heat with library sizes from lib_s2
  ratios_heat <- cbind(ratios_heat, Size = "")
  ratios_heat[, "Size"] <- lib_s2[match(rownames(ratios_heat), rownames(lib_s2)), "Size"]
  ratios_heat <- as.data.frame(ratios_heat)

    # Convert values of ratios_heat to numeric
  cols_to_convert <- c("total_POS", "total_GDNA", "total_gens", "total_NC", "total_ERCC",
                       "pos/gens", "gdna/gens", "nc/gens", "ERCC/gens", "median", "Size")
  ratios_heat[cols_to_convert] <- lapply(ratios_heat[cols_to_convert], as.numeric)
  str(ratios_heat)

    assign_01_QC <- function(valor, threshold) {
    ifelse(valor < threshold, 0, 1)
  }

  # Function to assign 0 or 1 according to library size value
  assign_01_size <- function(valor, threshold) {
    ifelse(valor > threshold, 0, 1)
  }

    # Create binary matrix for the heatmap
  bin_matrix <- matrix(0, nrow = nrow(ratios_heat), ncol = 6)
  for (i in 1:nrow(ratios_heat)) {
    bin_matrix[i, 1] <- assign_01_QC(ratios_heat[i, "pos/gens"], threshold_line_pos)
    bin_matrix[i, 2] <- assign_01_size(ratios_heat[i, "Size"], threshold_lib)
    bin_matrix[i, 3] <- assign_01_QC(ratios_heat[i, "nc/gens"], threshold_line_nc)
    bin_matrix[i, 4] <- assign_01_QC(ratios_heat[i, "gdna/gens"], threshold_line_gdna)
    bin_matrix[i, 5] <- assign_01_QC(ratios_heat[i, "ERCC/gens"], threshold_line_ercc)
    bin_matrix[i, 6] <- assign_01_size(ratios_heat[i, "median"], threshold_line_median)
  }

  # Row and column names
  rownames(bin_matrix) <- rownames(ratios_heat)
  colnames(bin_matrix) <- c("QC0", "QC1", "QC2", "QC3", "QC4", "QC5")

  # Convert the matrix to a data frame for ggplot2
  bin_df <- as.data.frame(bin_matrix)
  bin_df$Sample <- rownames(bin_df)
  # Melt the data frame
  bin_df_melted <- reshape2::melt(bin_df, id.vars = "Sample")

  # VIOLIN PLOT
  # Convert raw counts to TPM
    count2tpm <- function(countMat, idType = "Ensembl", org = "hsa",  source = "local", effLength = NULL, id = "id", gene_symbol = "symbol", length = "eff_length", check_data = FALSE){


  if(!org%in%c("hsa", "mmus")) stop(">>>== `org` must be hsa or mmus...")
  # requireNamespace("biomaRt")
  if(!is.matrix(countMat)){
    countMat<-as.matrix(countMat)
    countMat<-matrix(as.numeric(countMat), dim(countMat), dimnames = dimnames(countMat))
  }

  if(sum(is.na(countMat))>0|check_data){
    message(paste0("There are ", sum(is.na(countMat)) ," missing value in count matrix, these genes will be removed."))
    feas<-feature_manipulation(data = countMat, feature = rownames(countMat), is_matrix = T)
    countMat<-countMat[rownames(countMat)%in%feas,]

  }

  if(is.null(effLength) & source == "biomart"){
    datasets = paste0(c("hsapiens", "mmusculus", "btaurus", "cfamiliaris",
                        "ptroglodytes", "rnorvegicus", "sscrofa"), "_gene_ensembl")
    type = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol", "start_position", "end_position")
    if(org =="mmu") type[3] = "mgi_symbol"
    # listEnsemblArchives()
    # listMarts()
    # listAttributes()
    ds <- datasets[grepl(org, datasets)]

    ##
    ### In case it is not working you can find for a mirror place here: https://www.ensembl.org/info/about/mirrors.html
    mart <- biomaRt::useMart(host = "https://useast.ensembl.org", biomart = 'ENSEMBL_MART_ENSEMBL', dataset = ds)
    ensembl <- biomaRt::getBM(attributes=type, mart = mart)
    #######################################

    ensembl$Length <- abs(ensembl$end_position - ensembl$start_position)

    message(">>>--- This function is being optimised and we strongly recommend that you should set `source` as `local`....")
    #######################################
    if(toupper(idType) == "ENSEMBL"){

      len <- ensembl[match(rownames(countMat),ensembl$ensembl_gene_id), "Length"]
      rownames(countMat) = ensembl[match(rownames(countMat),ensembl$ensembl_gene_id), 3]
    }else if(toupper(idType) == "SYMBOL"){
      len <- ensembl[match(rownames(countMat), ensembl[,3]), "Length"]
    }else if(toupper(idType) == "ENTREZ"){
      len <- ensembl[match(rownames(countMat), ensembl[,2]), "Length"]
    }else{
      stop("Please input right type of gene name, such as `ensembl`, `entrez`, or `symbol` ...")
    }
  }


  if(source == "local" & tolower(idType) == "ensembl" & org == "hsa") {

    rownames(countMat) <- substring(rownames(countMat), 1, 15)
    data("anno_grch38", package = "IOBR")
    message(">>>--- Using variables (anno_grch38) and gene lengths (eff_length)  built into the IOBR package to perform TPM transformation")
    message(">>>--- The gene lengths (eff_length) was estimated by function `getGeneLengthAndGCContent` from EDASeq package with default parameters at 2023-02-10")

    length_ensembl<-anno_grch38[,c("id", "eff_length", "symbol")]
    length_ensembl<-length_ensembl[order(length_ensembl$eff_length, decreasing = T), ]

    countMat<- countMat[rownames(countMat)%in%length_ensembl$id,]

    if(dim(countMat)[1]==0) stop("Identifier of matrix is not match to references.")
    length_ensembl<-length_ensembl[length_ensembl$id%in%rownames(countMat),]
    len<- length_ensembl[match(rownames(countMat), length_ensembl$id), "eff_length"]
    rownames(countMat)<- length_ensembl[match(rownames(countMat),length_ensembl$id), 3]

    countMat <- matrix(as.numeric(countMat), dim(countMat), dimnames = dimnames(countMat))
  }else if(source == "local" & tolower(idType) == "entrez"  & org == "hsa"){


    message(">>>--- This is a fuzzy calculation. We recommend that users provide expression matrices with ENSEMBL as row names")
    message(">>>--- Using variables (anno_grch38) and gene lengths (eff_length)  built into the IOBR package to perform TPM transformation")
    message(">>>--- The gene lengths (eff_length) was estimated by function `getGeneLengthAndGCContent` from EDASeq package with default parameters at 2023-02-10")

    length_ensembl<-anno_grch38[,c("entrez", "eff_length", "symbol")]
    length_ensembl<-length_ensembl[order(length_ensembl$eff_length, decreasing = T), ]
    colnames(length_ensembl)[1]<-"id"
    length_ensembl<-length_ensembl[!duplicated(length_ensembl$id), ]

    countMat<-countMat[rownames(countMat)%in%length_ensembl$id, ]
    if(dim(countMat)[1]==0) stop("Identifier of matrix is not match to references.")
    length_ensembl<-length_ensembl[length_ensembl$id%in%rownames(countMat), ]
    len<- length_ensembl[match(rownames(countMat), length_ensembl$id), "eff_length"]
    rownames(countMat) <- length_ensembl[match(rownames(countMat),length_ensembl$id), 3]
    countMat<-matrix(as.numeric(countMat), dim(countMat), dimnames = dimnames(countMat))

  }else if(source == "local" & tolower(idType) == "symbol"  & org == "hsa"){

    message(">>>--- This is a fuzzy calculation. We recommend that users provide expression matrices with ENSEMBL as row names")
    message(">>>--- Using variables (anno_grch38) and gene lengths (eff_length)  built into the IOBR package to perform TPM transformation")
    message(">>>--- The gene lengths (eff_length) was estimated by function `getGeneLengthAndGCContent` from EDASeq package with default parameters at 2023-02-10")
    length_ensembl<-anno_grch38[, c("symbol", "eff_length", "gc")]
    length_ensembl<-length_ensembl[order(length_ensembl$eff_length, decreasing = T), ]

    colnames(length_ensembl)[1] <-"id"

    # print(head(length_ensembl))
    # print(head(countMat))

    length_ensembl <- length_ensembl[!duplicated(length_ensembl$id), ]
    countMat <- countMat[rownames(countMat)%in%length_ensembl$id, ]

    if(dim(countMat)[1]==0) stop("Identifier of matrix is not match to references.")
    length_ensembl<-length_ensembl[length_ensembl$id%in%rownames(countMat),]
    len<- length_ensembl[match(rownames(countMat), length_ensembl$id), "eff_length"]

    rownames(countMat)<- length_ensembl[match(rownames(countMat),length_ensembl$id), 1]
    countMat<-matrix(as.numeric(countMat), dim(countMat), dimnames = dimnames(countMat))
  }

  #######################################################################
  if(source == "local" & tolower(idType) == "ensembl" & org == "mmus") {

    message(">>>--- Using variables (anno_gc_vm32) and gene lengths (eff_length)  built into the IOBR package to perform TPM transformation")
    message(">>>--- The gene lengths (eff_length) was estimated by function `getGeneLengthAndGCContent` from EDASeq package with default parameters at 2023-02-10")

    length_ensembl<-anno_gc_vm32[,c("id", "eff_length", "symbol")]
    length_ensembl<-length_ensembl[order(length_ensembl$eff_length, decreasing = T), ]

    countMat<-countMat[rownames(countMat)%in%length_ensembl$id,]
    if(dim(countMat)[1]==0) stop("Identifier of matrix is not match to references.")
    length_ensembl<-length_ensembl[length_ensembl$id%in%rownames(countMat),]
    len<- length_ensembl[match(rownames(countMat), length_ensembl$id), "eff_length"]
    rownames(countMat)<- length_ensembl[match(rownames(countMat),length_ensembl$id), 3]
    countMat<-matrix(as.numeric(countMat), dim(countMat), dimnames = dimnames(countMat))
  }else if(source == "local" & tolower(idType) == "mgi"  & org == "mmus"){

    message(">>>--- This is a fuzzy calculation. We recommend that users provide expression matrices with ENSEMBL ID as row names")
    message(">>>--- Using variables (anno_gc_vm32) and gene lengths (eff_length)  built into the IOBR package to perform TPM transformation")
    message(">>>--- The gene lengths (eff_length) was estimated by function `getGeneLengthAndGCContent` from EDASeq package with default parameters at 2023-02-10")

    length_ensembl<-anno_gc_vm32[,c("mgi_id", "eff_length", "symbol")]
    length_ensembl<-length_ensembl[order(length_ensembl$eff_length, decreasing = T), ]
    colnames(length_ensembl)[1]<-"id"
    length_ensembl<-length_ensembl[!duplicated(length_ensembl$id), ]

    countMat<-countMat[rownames(countMat)%in%length_ensembl$id,]
    if(dim(countMat)[1]==0) stop("Identifier of matrix is not match to references.")
    length_ensembl<-length_ensembl[length_ensembl$id%in%rownames(countMat),]
    len<- length_ensembl[match(rownames(countMat), length_ensembl$id), "eff_length"]
    rownames(countMat) <- length_ensembl[match(rownames(countMat),length_ensembl$id), 3]
    countMat<-matrix(as.numeric(countMat), dim(countMat), dimnames = dimnames(countMat))

  }else if(source == "local" & tolower(idType) == "symbol"  & org == "mmus"){

    message(">>>--- This is a fuzzy calculation. We recommend that users provide expression matrices with ENSEMBL ID as row names")
    message(">>>--- Using variables (anno_gc_vm32) and gene lengths (eff_length)  built into the IOBR package to perform TPM transformation")
    message(">>>--- The gene lengths (eff_length) was estimated by function `getGeneLengthAndGCContent` from EDASeq package with default parameters at 2023-02-10")
    length_ensembl<-anno_gc_vm32[,c("symbol", "eff_length", "gc")]
    length_ensembl<-length_ensembl[order(length_ensembl$eff_length, decreasing = T), ]

    colnames(length_ensembl)[1]<-"id"
    length_ensembl<-length_ensembl[!duplicated(length_ensembl$id), ]

    countMat<-countMat[rownames(countMat)%in%length_ensembl$id,]
    if(dim(countMat)[1]==0) stop("Identifier of matrix is not match to references.")
    length_ensembl<-length_ensembl[length_ensembl$id%in%rownames(countMat),]
    len<- length_ensembl[match(rownames(countMat), length_ensembl$id), "eff_length"]
    rownames(countMat)<- length_ensembl[match(rownames(countMat),length_ensembl$id), 1]
    countMat<-matrix(as.numeric(countMat), dim(countMat), dimnames = dimnames(countMat))
  }

  #########################################################################
  if(!is.null(effLength)){
    effLength<-as.data.frame(effLength)
    colnames(effLength)[which(colnames(effLength)==id)]<-"id"
    colnames(effLength)[which(colnames(effLength)==length)]<-"eff_length"
    effLength<-effLength[!duplicated(effLength$id),]

    countMat<-as.matrix(countMat)
    countMat<-countMat[rownames(countMat)%in%effLength$id, ]
    effLength<-effLength[effLength$id%in%rownames(countMat), ]

    if(id!= gene_symbol){
      # countMat<-as.matrix(countMat)
      colnames(effLength)[which(colnames(effLength)==gene_symbol)]<-"gene_symbol"
      rownames(countMat)<- effLength[match(rownames(countMat),effLength$id), "gene_symbol"]

    }else{
      # countMat<-as.matrix(countMat)
      effLength$gene_symbol<-effLength$id
      # colnames(effLength)[which(colnames(effLength)==gene_symbol)]<-"gene_symbol"
      rownames(countMat)<- effLength[match(rownames(countMat),effLength$id), "gene_symbol"]
    }

    len<- effLength[match(rownames(countMat), effLength[,"gene_symbol"]), "eff_length"]

  }

  na_idx <- which(is.na(len))
  if(length(na_idx)>0){
    warning(paste0(">>>--- Omit ", length(na_idx), " genes of which length is not available !"))
    countMat <- countMat[!is.na(len), ]
    len = len[!is.na(len)]
  }
  #####################################
  tmp <- countMat / c(len/1000) # (`per million` scaling factor)
  TPM <- 1e6 * t(t(tmp) / colSums(tmp))
  TPM <- TPM[!is.na(rownames(TPM)),]
  TPM <- TPM[!rownames(TPM)==" ",]

  # TPM <- rownames_to_column(as.data.frame(TPM), var = "symbol")
  symbol.id = rownames(TPM)
  TPM = as.data.frame(TPM)
  TPM$symbol = symbol.id

  TPM <- remove_duplicate_genes(eset = TPM, column_of_symbol = "symbol")
  # TPM <- TPM[,!is.na(colnames(TPM))]
  # TPM <- TPM[,!colnames(TPM)==" "]
  return(TPM)
}
  tpm_counts <- count2tpm(counts_data,
                          idType = "Symbol",
                          org = "hsa",
                          source = "biomart")

  # TPM data formatting
  tpm_counts$Gene <- rownames(tpm_counts)
  tpm_long <- reshape2::melt(tpm_counts, id.vars = "Gene", variable.name = "Sample", value.name = "TPM")
  tpm_long$LogTPM <- log1p(tpm_long$TPM)

  # Calculate the 95th percentile threshold for each sample
  #percentile_95 <- aggregate(LogTPM ~ Sample, data = tpm_long, FUN = function(x) quantile(x, 0.95))
  percentile_95 <- dplyr::summarise(dplyr::group_by(tpm_long, Sample),
    percentile_95 = quantile(LogTPM, 0.95))


  colnames(percentile_95)[2] <- "Threshold"
  tpm_with_threshold <- merge(tpm_long, percentile_95, by = "Sample")

  # Filter data based on the 95th percentile threshold
  tpm_filtered <- tpm_with_threshold[tpm_with_threshold$LogTPM < tpm_with_threshold$Threshold, ]

  # Function to create violin plot
  create_violin_plot <- function(data, title) {
    ggplot2::ggplot(data, ggplot2::aes(x = Sample, y = LogTPM)) +
      ggplot2::geom_violin(trim = FALSE, fill = "#4793AF", color = "black") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
                     panel.background = ggplot2::element_rect(fill = "white"),
                     plot.background = ggplot2::element_rect(fill = "white"),
                     panel.grid.major = ggplot2::element_line(color = "gray80"),
                     panel.grid.minor = ggplot2::element_line(color = "gray90")) +
      ggplot2::labs(title = title,
                    x = "Sample",
                    y = "Log-Transformed TPM")
  }

  # Check if there are more than 50 unique samples and split if needed
  if (length(unique(tpm_filtered$Sample)) > 50) {
    samples <- unique(tpm_filtered$Sample)
    half <- ceiling(length(samples) / 2)
    subset1 <- samples[1:half]
    subset2 <- samples[(half + 1):length(samples)]
    data1 <- tpm_filtered[tpm_filtered$Sample %in% subset1, ]
    data2 <- tpm_filtered[tpm_filtered$Sample %in% subset2, ]
    p4 <- create_violin_plot(data1, "Distribution of Log-Transformed TPM (Up to 95th Percentile) - Part 1")
    p5 <- create_violin_plot(data2, "Distribution of Log-Transformed TPM (Up to 95th Percentile) - Part 2")
    combined_plot <- ggpubr::ggarrange(p4, p5, ncol = 1, nrow = 2)
  } else {
    p6 <- create_violin_plot(tpm_filtered, "Distribution of Log-Transformed TPM (Up to 95th Percentile)")
    print(summary(tpm_filtered))
    combined_plot <- p6
  }

  # Save the violin plot to a PDF
  pdf("QC_plots_violin_plot.pdf", width = 14, height = 10)
  print(combined_plot)
  dev.off()


  cat("\033[32mViolin plot saved as 'QC_plots_violin_plot.pdf'\033[0m\n")


  pdf("QC_plots.pdf")

  # Positive controls
  max_value <- max(ratios$`pos/gens`, threshold_line_pos)
  min_size <- min(ratios$`pos/gens`, threshold_line_pos)

  colores_pos <- ifelse(ratios$`pos/gens` <= threshold_line_pos, "#4793AF",
                        ifelse(ratios$`pos/gens` <= threshold_superior_pos, "#FFC470", "red"))
  plot(ratios$`pos/gens`, xlab = "", ylab = "pos/gens", col = colores_pos,
       xaxt = "n", pch = 19, main = "Positive control 4% (QC0)", ylim = c(0, max_value))
  axis(1, at = 1:nrow(ratios), labels = rownames(ratios), las = 2, cex.axis = 0.8)
  abline(h = threshold_line_pos, col = "red")


  # Library size
  max_size <- max(lib_s2$Size, threshold_lib)
  min_size <- min(lib_s2$Size, threshold_lib)
    colores <- ifelse(lib_s2$Size < threshold_inferior_lib, "red",
                    ifelse(lib_s2$Size <= threshold_lib, "#FFC470", "#4793AF"))
  plot(lib_s2$Size, xlab = "", ylab = "Library Size", col = colores,
       xaxt = "n", pch = 19, main = "Library Size per Sample (QC1)", cex.axis = 0.8,
       ylim = c(min_size, max_size))
  axis(1, at = 1:length(lib_s2$Sample), labels = lib_s2$Sample, las = 2, cex.axis = 0.8)
  abline(h = threshold_lib, col = "red")


  # Negative controls
  max_value <- max(ratios$`nc/gens`, threshold_line_nc)
  min_value <- min(ratios$`nc/gens`, threshold_line_nc)
  colores_nc <- ifelse(ratios$`nc/gens` <= threshold_line_nc, "#4793AF",
                       ifelse(ratios$`nc/gens` <= threshold_superior_nc, "#FFC470", "red"))
  plot(ratios$`nc/gens`, xlab = "", ylab = "nc/gens", col = colores_nc,
       xaxt = "n", pch = 19, main = "Negative Control (QC2)", ylim = c(min_value, max_value))
  axis(1, at = 1:nrow(ratios), labels = rownames(ratios), las = 2, cex.axis = 0.8)
  abline(h = threshold_line_nc, col = "red")


  # Genomic DNA
  max_value <- max(ratios$`gdna/gens`, threshold_line_gdna)
  min_size <- min(ratios$`gdna/gens`, threshold_line_gdna)
  colores_gdna <- ifelse(ratios$`gdna/gens` <= threshold_line_gdna, "#4793AF",
                         ifelse(ratios$`gdna/gens` <= threshold_superior_gdna, "#FFC470", "red"))
  plot(ratios$`gdna/gens`, xlab = "", ylab = "gdna/gens", col = colores_gdna,
       xaxt = "n", pch = 19, main = "Genomic DNA (QC3)", ylim = c(0, max_value))
  axis(1, at = 1:nrow(ratios), labels = rownames(ratios), las = 2, cex.axis = 0.8)
  abline(h = threshold_line_gdna, col = "red")

  # ERCC
  max_value <- max(ratios$`ERCC/gens`, threshold_line_ercc)
  min_size <- min(ratios$`ERCC/gens`, threshold_line_ercc)

  colores_ercc <- ifelse(ratios$`ERCC/gens` <= threshold_line_ercc, "#4793AF",
                         ifelse(ratios$`ERCC/gens` <= threshold_superior_ercc, "#FFC470", "red"))
  plot(ratios$`ERCC/gens`, xlab = "", ylab = "ERCC/gens", col = colores_ercc,
       xaxt = "n", pch = 19, main = "ERCC (QC4)", ylim = c(0, max_value))
  axis(1, at = 1:nrow(ratios), labels = rownames(ratios), las = 2, cex.axis = 0.8)
  abline(h = threshold_line_ercc, col = "red")


  ## Median
  max_value <- max(ratios$median, threshold_line_median)
  min_size <- min(ratios$median, threshold_line_median)
  colores_med <- ifelse(ratios$median < threshold_inferior_median, "red",
                        ifelse(ratios$median <= threshold_line_median, "#FFC470", "#4793AF"))
  plot(ratios$median, xlab = "", ylab = "Median", col = colores_med,
       xaxt = "n", pch = 19, main = "Median (QC5)", ylim = c(0, max_value))
  axis(1, at = 1:nrow(ratios), labels = rownames(ratios), las = 2, cex.axis = 0.8)
  abline(h = threshold_line_median, col = "red")
    dev.off()

    cat("\033[32mQC plots saved as 'plot_QC.pdf'\033[0m\n")

  # Create the heatmap with ggplot2
  a<- ggplot2::ggplot(bin_df_melted, ggplot2::aes(x = variable, y = Sample, fill = factor(value))) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_manual(values = c("0" = "#FFF9D0", "1" = "red"), labels = c("OK", "Outlier")) +
    ggplot2::labs(x = "QC Metrics", y = " ", fill = "QC Status") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
          axis.text.y = ggplot2::element_text(size = 7),
          legend.position = "bottom")

  pdf("HTG_QC_heatmap.pdf", width = 10, height = 14)
  print(a)
      dev.off()

  cat("\033[32mHeatmap PLOTS SAVED AS 'HTG_QC_heatmap.pdf'\033[0m\n")

  rows_with_1 <- suppressWarnings(rownames(bin_matrix)[apply(bin_matrix, 1, any)])
  cat("\033[32m                              ***\033[0m\n")
  cat(paste("\033[32mThese are the samples plotted at least once in the heatmap:  \033[0m\n"))
  cat("The number of samples that are outliers is:", length(rows_with_1), "\n")

  cat(rows_with_1)
  return(rows_with_1)
}


