#' HTG_TME
#'
#' @description This function carries out a comprehensive TME analysis involving TPM normalization and TME using multiple methods (EPIC, quanTIseq, and xCell). It produces several output files containing the results of the normalization and TME processes.

#' @param outliers Vector of outlier sample IDs to be removed from analysis.
#' @param pattern Regular expression pattern to identify outlier sample IDs. For HTG, this could be "^NC-|^POS-|^GDNA-|^ERCC-". If NULL, the pattern will not be applied.
#' @param counts_data Count data matrix of gene expression.
#' @param AnnotData DataFrame containing sample metadata.
#' @param design_formula Formula specifying the design of the analysis.
#' @param DEA Preprocessed DESeqDataSet for differential expression analysis.
#' @param generate_volcano Logical; whether to generate volcano plots.
#' @param remove_outliers Logical; whether to remove outlier samples.
#'
#' @return The function outputs multiple CSV files (TPM normalized counts and Results of TME methods EPIC, quanTIseq, and xCell), PDF with graphics and return a dataframe with the three TME methods.
#'
#' @export
#'
#' @import immunedeconv
#'
#' @examples
#' TME_results <- HTG_TME(
#'   outliers = outliers_tutorial,
#'   pattern = "^NC-|^POS-|^GDNA-|^ERCC-",
#'   counts_data = counts_data_tutorial,
#'   AnnotData = AnnotData_tutorial,
#'   design_formula = "HPV_status",
#'   remove_outliers = TRUE,
#'   DEA = NULL
#' )
#'
#' @name HTG_TME
#'
#'
#'
utils::globalVariables(c(".data", "mean_value", "Sample", "Cell_Type", "Value", "Average", "Fraction", "shapiro_test"))

HTG_TME <- function(outliers, pattern = NULL, counts_data, AnnotData, design_formula = NULL ,
                            DEA = NULL, generate_volcano = TRUE, remove_outliers = TRUE) {

  if (!requireNamespace("IOBR", quietly = TRUE)) {
    stop("Package 'IOBR' is required but not installed.")
  }
  if (!requireNamespace("immunedeconv", quietly = TRUE)) {
    stop("Package 'immunedeconv' is required but not installed.")
  }

  if (!requireNamespace("EPIC", quietly = TRUE)) {
    stop("Package 'EPIC' is required but not installed.")
  }

  if (!requireNamespace("xCell", quietly = TRUE)) {
    stop("Package 'xCell' is required but not installed.")
  }
  utils::data("xCell.data", package = "xCell", envir = environment())


  if (!is.null(pattern)) {
    # Remove outliers based on pattern
    filtered <- subset(counts_data, !grepl(pattern, rownames(counts_data)))
  } else {
    # If no pattern is provided, do not filter based on pattern
    filtered <- counts_data
  }

  if (remove_outliers) {
    # Remove columns corresponding to outliers
    counts_data <- filtered[, !colnames(filtered) %in% outliers]
    AnnotData <- AnnotData[!AnnotData[["id"]] %in% outliers, ]
  } else {
    counts_data <- counts_data
    AnnotData <- AnnotData
  }

  print(paste0("Inicial gene number: ", dim(counts_data)[1]))
  ## Normalización TPM
  cat("\033[32mTPM normalization performed and stored on tpm_counts.csv\033[0m\n")
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
    mart <- biomaRt::useMart(host = "https://www.ensembl.org", biomart = 'ENSEMBL_MART_ENSEMBL', dataset = ds)
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
  suppressWarnings({
    tpm_counts <- count2tpm(counts_data,
                            idType = "Symbol",
                            org = "hsa",
                            source = "biomart")
    write.csv(tpm_counts, "tpm_counts.csv")
  })

  # Se almacenan los genes omitidos en la normalización TPM
  genes_omitidos <- base::setdiff(rownames(counts_data), rownames(tpm_counts))
  print(paste0("Number of genes omitted during TPM normalization due to their length not being available in Biomart:   ", dim(counts_data)[1]-dim(tpm_counts)[1]))
  tpm_counts <- as.data.frame(tpm_counts)

  if (!is.null(DEA)) {
    cat("\033[32mWe are going to use information from DEA\033[0m\n")
    dds <- DESeq2::estimateSizeFactors(dds)
    normalized_counts <- DESeq2::counts(dds, normalized = TRUE)
  } else {
    cat("Performing normalization.\n")
    design_formul <- as.formula("~ 1")
    colnames(AnnotData) <- gsub(" ", "_", colnames(AnnotData))
    AnnotData <- AnnotData[order(AnnotData$id), ]
    counts_data <- counts_data[, order(colnames(counts_data))]
    if (!identical(colnames(counts_data), AnnotData$id)) {
      stop("Column names of counts_data and IDs in AnnotData do not match.")
    }
  }

  ## Deconvolución
  imm_epic <- immunedeconv::deconvolute(tpm_counts, method = "epic")
  imm_qti <- immunedeconv::deconvolute(tpm_counts, method = "quantiseq")
  imm_xcell <- immunedeconv::deconvolute(tpm_counts, method = "xcell")

  cat("\033[32mresults of the devonvolution will be stored in imm_epic.csv, imm_qti.csv and imm_xcell.csv \033[0m\n")
  write.csv(imm_epic, file = "imm_epic.csv")
  write.csv(imm_qti, file = "imm_qti.csv")
  write.csv(imm_xcell, file = "imm_xcell.csv")

  # Transponer y cambiar colnames
  std.im.df <- function(imm_df){
    imm_df <- as.data.frame(t(imm_df))
    celltype_imm <- imm_df[1,]
    imm_df <- imm_df[-1,]
    colnames(imm_df) <- celltype_imm
    imm_df[-1,]
    rn_imm <- rownames(imm_df)
    imm_df <- as.data.frame(sapply(imm_df, as.numeric))
    rownames(imm_df) <- rn_imm
    return(imm_df)
  }

  imm_epic <- std.im.df(imm_epic)
  imm_qti <- std.im.df(imm_qti)
  imm_xcell <- std.im.df(imm_xcell)

  # Se incluye la variable AnnotData$design_formula en los dataframes
  cat("\033[32mAre they in the same order?\033[0m\n")

  rownames(AnnotData)<- AnnotData$id
  AnnotData <- AnnotData[order(rownames(AnnotData)), ]
  imm_epic <- imm_epic[order(rownames(imm_epic)), ]
  imm_qti <- imm_qti[order(rownames(imm_qti)), ]
  imm_xcell <- imm_xcell[order(rownames(imm_xcell)), ]
  AnnotData$id
  rownames(imm_epic)
  rownames(imm_qti)
  rownames(imm_xcell)

  cat("\033[32mHave to be true.\033[0m\n")
  cat("\033[32mEPIC\033[0m\n")
  print(all(rownames(AnnotData)==rownames(imm_epic)))
  cat("\033[32mqti\033[0m\n")
  print(all(rownames(AnnotData)==rownames(imm_qti)))
  cat("\033[32mxcell\033[0m\n")
  print(all(rownames(AnnotData)==rownames(imm_xcell)))

  imm_epic[[design_formula]] <- factor(AnnotData[[design_formula]])
  imm_qti[[design_formula]] <- factor(AnnotData[[design_formula]])
  imm_xcell[[design_formula]] <- factor(AnnotData[[design_formula]])

  replace_space_with_underscore <- function(df) {
    colnames(df) <- gsub(" ", "_", colnames(df))
    return(df)
  }
  imm_epic <- replace_space_with_underscore(imm_epic)
  imm_qti <- replace_space_with_underscore(imm_qti)
  imm_xcell <- replace_space_with_underscore(imm_xcell)

  ############# imm_epic
  design_formula_sym <- rlang::sym(design_formula)
  # Vector with column names except the last one
  column_names <- names(imm_epic)[-ncol(imm_epic)]
  # Loop to generate plots for each column

  pdf("plots_imm_EPIC.pdf")
  for (col_name in column_names) {
    # Calculate means by design formula

    grouped_data <- dplyr::group_by(imm_epic, !!design_formula_sym)
    means_df <- dplyr::summarize(grouped_data, mean_value = mean(.data[[col_name]] * 100, na.rm = TRUE), .groups = "drop")
    # Perform t-test or ANOVA
    formula <- as.formula(paste0("`", col_name, "` ~ ", rlang::as_label(design_formula_sym)))
    p_value <- tryCatch({
      if (dplyr::n_distinct(imm_epic[[rlang::as_label(design_formula_sym)]]) == 2) {
        t.test(formula, data = imm_epic)$p.value
      } else {
        summary(stats::aov(formula, data = imm_epic))[[1]]$`Pr(>F)`[1]
      }
    }, error = function(e) {
      NA  # In case of error, return NA for p-value
    })
    # Generate plot
    plot <- ggplot2::ggplot(imm_epic, ggplot2::aes(x = !!design_formula_sym , y = .data[[col_name]] * 100,
                                                   fill = !!design_formula_sym, color = !!design_formula_sym)) +
      ggplot2::geom_jitter(alpha = 1, width = 0.3, height = 0) +
      ggplot2::geom_boxplot(fill = "white", alpha = 0.5, outlier.alpha = 1) +
      ggplot2::geom_point(data = means_df, ggplot2::aes(x = !!design_formula_sym, y = mean_value, fill = !!design_formula_sym),
                          shape = 22, color = "black", size = 3, stroke = 1.5,
                          show.legend = F) +
      ggplot2::labs(x = NULL, y = "Abundance (%)", title = col_name, subtitle = "imm_epic") +
      ggplot2::scale_y_continuous(labels = scales::percent_format(scale = 1),
                                  limits = c(0, NA)) +
      ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
      ggplot2::annotate("text", x = Inf, y = Inf, label = paste("p-value:", format(p_value, digits = 3)),
                        hjust = 1.1, vjust = 1.1, size = 5, color = "red")
    print(plot)
  }
  dev.off()

  ############# imm_qti
  # Vector with column names except the last one
  column_names <- names(imm_qti)[-ncol(imm_qti)]

  pdf("plots_imm_qti.pdf")
  for (col_name in column_names) {
    # Calculate means by design formula
    grouped_data <- dplyr::group_by(imm_qti, !!design_formula_sym)
    means_df <- dplyr::summarize(grouped_data, mean_value = mean(.data[[col_name]] * 100, na.rm = TRUE), .groups = "drop")

    # Perform t-test or ANOVA
    formula <- as.formula(paste0("`", col_name, "` ~ ", rlang::as_label(design_formula_sym)))
    p_value <- tryCatch({
      if (dplyr::n_distinct(imm_qti[[rlang::as_label(design_formula_sym)]]) == 2) {
        t.test(formula, data = imm_qti)$p.value
      } else {
        summary(stats::aov(formula, data = imm_qti))[[1]]$`Pr(>F)`[1]
      }
    }, error = function(e) {
      NA  # In case of error, return NA for p-value
    })

    # Generate plot
    plot <- ggplot2::ggplot(imm_qti, ggplot2::aes(x = !!design_formula_sym , y = .data[[col_name]] * 100,
                                                  fill = !!design_formula_sym, color = !!design_formula_sym)) +
      ggplot2::geom_jitter(alpha = 1, width = 0.3, height = 0) +
      ggplot2::geom_boxplot(fill = "white", alpha = 0.5, outlier.alpha = 1) +
      ggplot2::geom_point(data = means_df, ggplot2::aes(x = !!design_formula_sym, y = mean_value, fill = !!design_formula_sym),
                          shape = 22, color = "black", size = 3, stroke = 1.5,
                          show.legend = F) +
      ggplot2::labs(x = NULL, y = "Abundance (%)", title = col_name, subtitle = "imm_qti") +
      ggplot2::scale_y_continuous(labels = scales::percent_format(scale = 1),
                                  limits = c(0, NA)) +
      ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
      ggplot2::annotate("text", x = Inf, y = Inf, label = paste("p-value:", format(p_value, digits = 3)),
                        hjust = 1.1, vjust = 1.1, size = 5, color = "red")
    print(plot)
  }
  dev.off()

  ############# imm_xcell
  # Vector with column names except the last one
  column_names <- names(imm_xcell)[-ncol(imm_xcell)]
  pdf("plots_imm_xcell.pdf")
  for (col_name in column_names) {
    # Calculate means by design formula
    grouped_data <- dplyr::group_by(imm_xcell, !!design_formula_sym)
    means_df <- dplyr::summarize(grouped_data, mean_value = mean(.data[[col_name]] * 100, na.rm = TRUE), .groups = "drop")

    # Perform t-test or ANOVA
    formula <- as.formula(paste0("`", col_name, "` ~ ", rlang::as_label(design_formula_sym)))
    p_value <- tryCatch({
      if (dplyr::n_distinct(imm_xcell[[rlang::as_label(design_formula_sym)]]) == 2) {
        t.test(formula, data = imm_xcell)$p.value
      } else {
        summary(stats::aov(formula, data = imm_xcell))[[1]]$`Pr(>F)`[1]
      }
    }, error = function(e) {
      NA  # In case of error, return NA for p-value
    })
    # Generate plot
    plot <- ggplot2::ggplot(imm_xcell, ggplot2::aes(x = !!design_formula_sym , y = .data[[col_name]] * 100,
                                                    fill = !!design_formula_sym, color = !!design_formula_sym)) +
      ggplot2::geom_jitter(alpha = 1, width = 0.3, height = 0) +
      ggplot2::geom_boxplot(fill = "white", alpha = 0.5, outlier.alpha = 1) +
      ggplot2::geom_point(data = means_df, ggplot2::aes(x = !!design_formula_sym, y = mean_value, fill = !!design_formula_sym),
                          shape = 22, color = "black", size = 3, stroke = 1.5,
                          show.legend = F) +
      ggplot2::labs(x = NULL, y = "Abundance (%)", title = col_name, subtitle = "imm_xcell") +
      ggplot2::scale_y_continuous(labels = scales::percent_format(scale = 1),
                                  limits = c(0, NA)) +
      ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
      ggplot2::annotate("text", x = Inf, y = Inf, label = paste("p-value:", format(p_value, digits = 3)),
                        hjust = 1.1, vjust = 1.1, size = 5, color = "red")
    print(plot)
  }
  dev.off()
  # Composición celular del TME i por grupo
  # Function to generate the bar plot for each sample
  plot_bar <- function(df, paleta, titulo, legend.position) {
    # Convertir las filas en una columna llamada "Sample"
    df <- tibble::rownames_to_column(df, var = "Sample")
    df <- tidyr::pivot_longer(
      df, cols = colnames(df)[2:(ncol(df) - 1)],
      names_to = "Cell_Type", values_to = "Value")
    df <- dplyr::mutate(df,
                        Sample = factor(Sample, levels = rev(unique(Sample))),
                        Cell_Type = factor(Cell_Type, levels = rev(unique(Cell_Type))))

    p <- ggplot2::ggplot(df, ggplot2::aes(x = Sample, y = Value, fill = Cell_Type)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::labs(title = titulo,
                    x = "Samples",
                    y = "Cell Fraction (%)") +
      ggplot2::coord_flip() +
      ggplot2::guides(fill = ggplot2::guide_legend(reverse = TRUE)) +
      ggplot2::scale_fill_manual(values = paleta) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = legend.position,
                     axis.text.y = ggplot2::element_text(size = 5)) +
      ggplot2::scale_y_continuous(labels = scales::percent)

    return(p)
  }

  plot_bar_group <- function(df, paleta, titulo, design_formula, legend_position = "right") {
    suppressWarnings({
      design_formula_sym <- rlang::sym(design_formula)
      niveles_tipo_cel <- colnames(df)[1:(ncol(df) - 1)]
      df_rownames <- tibble::rownames_to_column(df, var = "Sample")
      df_long <- tidyr::pivot_longer(
        df_rownames,
        cols = niveles_tipo_cel,
        names_to = "Cell_Type",
        values_to = "Value")
      df_grouped <- dplyr::group_by(df_long, !!design_formula_sym, Cell_Type)
      df_summarised <- dplyr::summarise(df_grouped,
                                        Average = mean(Value, na.rm = TRUE),
                                        .groups = "drop")
      df_ungrouped <- dplyr::ungroup(df_summarised)
      promedios <- dplyr::mutate(df_ungrouped,
                                 !!design_formula_sym := factor(!!design_formula_sym, levels = rev(unique(!!design_formula_sym))),
                                 Cell_Type = factor(Cell_Type, levels = rev(niveles_tipo_cel)))
      p <- ggplot2::ggplot(promedios, ggplot2::aes(x = !!design_formula_sym, y = Average, fill = Cell_Type)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::labs(title = titulo,
                      x = "Samples",
                      y = "Cell Fraction (%)") +
        ggplot2::coord_flip() +
        ggplot2::guides(fill = ggplot2::guide_legend(reverse = TRUE)) +
        ggplot2::scale_fill_manual(values = paleta) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = legend_position,
                       axis.text.y = ggplot2::element_text(size = 5)) +
        ggplot2::scale_y_continuous(labels = scales::percent)
      return(p)
    })
  }

  # Function to combine both plots into one
  plot_combined <- function(df, paleta, titulo_individual, titulo_grupo, design_formula, legend_position = "right") {
    p1 <- plot_bar(df, paleta, titulo_individual, legend_position)
    p2 <- plot_bar_group(df, paleta, titulo_grupo, design_formula, legend_position)

    combined_plot <- ggpubr::ggarrange(p1, p2, ncol = 1, nrow = 2, heights = c(1, 1))
    return(combined_plot)
  }

  # Example usage with extended palette for larger datasets
  paleta_imm <- c("grey95","#FB8072","#FFED6F","#6F6C87","#94DFD1","#FDB462", "#B3DE69", "#FFB1D9")
  paleta_qti <- c("grey95","#B3DE69","#6F6C87","#94DFD1","#FDB462", "#FB8072","#FFFFB3","#8BB07A","#FFED6F","#80B1D3","#FFB1D9")
  paleta_extendida <- c("grey95","#2CA02C","#E6F8E0","#7F7F7F","#FF8000","#DF7401","#6F6C80","#F5A9E1",
                        "#F7D358","#6F6C99","#FB8072","#FFFFB3","#F5A9F2",
                        "#8BB07A","#E377C2","#FE2E2E","#FFED6F","#1F77B4","#80B1D3","#0489B1","#FF0000","#8BB07A",
                        "#0000FF","#D6616B","#58FA82","#98DF8A",
                        "#8C564B", "#C49C94","#CEE3F6","#E0F8EC","#58FAF4", "#ADD8E6","#9EDAE5","#94DFD1","#FFA500",
                        "#FF7F0E","#FDB462","#FFB1D9","#B3DE69")
  # Assuming 'imm_epic', 'imm_qti', and 'imm_xcell' dataframes are already loaded
  combined_plot_EPIC <- plot_combined(imm_epic, paleta_imm, "EPIC Individual", "EPIC Average", design_formula, "right")
  combined_plot_quanTIseq <- plot_combined(imm_qti, paleta_qti, "quanTIseq Individual", "quanTIseq Average", design_formula, "right")
  # Mensaje de atención en rojo
  cat("\033[33mAttention: plot_cell_fraction_Average_cell_fraction will work for all methods. For xCell, it will display Enrichment Scores instead of percentages.\033[0m\n")

  plot_bar <- function(df, paleta, titulo, legend.position) {
    # Convertir las filas en una columna llamada "Sample"
    df <- tibble::rownames_to_column(df, var = "Sample")
    df <- tidyr::pivot_longer(
      df, cols = colnames(df)[2:(ncol(df) - 1)],
      names_to = "Cell_Type", values_to = "Value")
    df <- dplyr::mutate(df,
                        Sample = factor(Sample, levels = rev(unique(Sample))),
                        Cell_Type = factor(Cell_Type, levels = rev(unique(Cell_Type))))

    p <- ggplot2::ggplot(df, ggplot2::aes(x = Sample, y = Value, fill = Cell_Type)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::labs(title = titulo,
                    x = "Samples",
                    y = "Enrichment Scores") +
      ggplot2::coord_flip() +
      ggplot2::guides(fill = ggplot2::guide_legend(reverse = TRUE)) +
      ggplot2::scale_fill_manual(values = paleta) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = legend.position,
                     axis.text.y = ggplot2::element_text(size = 5)) +
      ggplot2::scale_y_continuous(labels = scales::percent)

    return(p)
  }

  plot_bar_group <- function(df, paleta, titulo, design_formula, legend_position = "right") {
    suppressWarnings({
      design_formula_sym <- rlang::sym(design_formula)
      niveles_tipo_cel <- colnames(df)[1:(ncol(df) - 1)]
      df_rownames <- tibble::rownames_to_column(df, var = "Sample")
      df_long <- tidyr::pivot_longer(
        df_rownames,
        cols = niveles_tipo_cel,
        names_to = "Cell_Type",
        values_to = "Value")
      df_grouped <- dplyr::group_by(df_long, !!design_formula_sym, Cell_Type)
      df_summarised <- dplyr::summarise(df_grouped,
                                        Average = mean(Value, na.rm = TRUE),
                                        .groups = "drop")
      df_ungrouped <- dplyr::ungroup(df_summarised)
      promedios <- dplyr::mutate(df_ungrouped,
                                 !!design_formula_sym := factor(!!design_formula_sym, levels = rev(unique(!!design_formula_sym))),
                                 Cell_Type = factor(Cell_Type, levels = rev(niveles_tipo_cel)))
      p <- ggplot2::ggplot(promedios, ggplot2::aes(x = !!design_formula_sym, y = Average, fill = Cell_Type)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::labs(title = titulo,
                      x = "Samples",
                      y = "Enrichment Scores") +
        ggplot2::coord_flip() +
        ggplot2::guides(fill = ggplot2::guide_legend(reverse = TRUE)) +
        ggplot2::scale_fill_manual(values = paleta) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = legend_position,
                       axis.text.y = ggplot2::element_text(size = 5)) +
        ggplot2::scale_y_continuous(labels = scales::percent)
      return(p)
    })
  }

  # Function to combine both plots into one
  plot_combined <- function(df, paleta, titulo_individual, titulo_grupo, design_formula, legend_position = "right") {
    p1 <- plot_bar(df, paleta, titulo_individual, legend_position)
    p2 <- plot_bar_group(df, paleta, titulo_grupo, design_formula, legend_position)

    combined_plot <- ggpubr::ggarrange(p1, p2, ncol = 1, nrow = 2, heights = c(1, 1))
    return(combined_plot)
  }
  combined_plot_xCell <- plot_combined(imm_xcell, paleta_extendida, "xCell Individual", "xCell Average", design_formula, "right")

  pdf("plot_cell_fraction_Average_cell_fraction_EPIC.pdf", width = 11, height = 14)
  print(combined_plot_EPIC)
  dev.off()

  pdf("plot_cell_fraction_Average_cell_fraction_quanTIseq.pdf", width = 11, height = 14)
  print(combined_plot_quanTIseq)
  dev.off()

  pdf("plot_cell_fraction_Average_cell_fraction_xCell.pdf", width = 11, height = 14)
  print(combined_plot_xCell)
  dev.off()

  ########################
  trans_formato_largo <- function(df, design_formula_sym) {
    df_largo <- tidyr::pivot_longer(df, cols = -dplyr::all_of(design_formula_sym), names_to = "Cell_Type",
                                    values_to = "Fraction")
    return(df_largo)
  }

  # Function to perform Shapiro-Wilk normality test
  prueba_norm <- function(df, design_formula_sym) {
    df_largo <- trans_formato_largo(df, design_formula_sym)

    df_grouped <- dplyr::group_by(df_largo, Cell_Type)
    df_summarised <- dplyr::summarise(df_grouped, shapiro_test = list(stats::shapiro.test(Fraction)))
    resultados_normalidad <- dplyr::mutate(df_summarised, p.value = purrr::map_dbl(shapiro_test, "p.value"))

    return(resultados_normalidad)
  }

  # Function to check if all values in a column are the same
  check_column_equal <- function(column) {
    all_equal <- length(unique(column)) == 1
    return(all_equal)
  }

  # Function to filter valid columns
  filter_valid_columns <- function(df, design_formula_sym, dataset_name) {
    cat("\033[32mChecking if All Values in Each Column Are Equal in dataset:", dataset_name, "\033[0m\n")

    equality_results <- base::sapply(df, check_column_equal)
    unequal_columns_count <- base::sum(!equality_results)

    if (unequal_columns_count > 0) {
      cat("\033[33mIn dataset", dataset_name, "the following columns have the same value in all rows and will be excluded from the analysis:\033[0m\n")
      print(names(df)[equality_results])
    } else {
      cat("\033[32mAll columns in dataset", dataset_name, "have variability. No columns to exclude.\033[0m\n")
    }

    valid_columns <- !equality_results
    df_filtered <- df[, valid_columns]

    return(df_filtered)
  }

  # Function to perform parametric tests (t-test and ANOVA)
  perform_parametric_tests <- function(df, design_formula_sym) {
    df_largo <- trans_formato_largo(df, design_formula_sym)
    cell_types <- unique(df_largo$Cell_Type)

    results <- list()

    for (cell_type in cell_types) {

      df_cell_type <- dplyr::filter(df_largo, Cell_Type == cell_type)
      groups <- unique(df_cell_type[[design_formula_sym]])


      if (length(groups) == 2) {
        # Perform t-test
        t_test_result <- t.test(Fraction ~ df_cell_type[[design_formula_sym]], data = df_cell_type)
        results[[cell_type]] <- list(Test = "t-test", p.value = t_test_result$p.value)
      } else {
        # Perform ANOVA
        anova_result <- stats::aov(Fraction ~ df_cell_type[[design_formula_sym]], data = df_cell_type)
        p_value <- summary(anova_result)[[1]]$`Pr(>F)`[1]
        results[[cell_type]] <- list(Test = "ANOVA", p.value = p_value)
      }
    }

    results_df <- dplyr::bind_rows(lapply(names(results), function(cell_type) {
      result <- results[[cell_type]]
      tibble::tibble(Cell_Type = cell_type, Test = result$Test, p.value = result$p.value)
    }))

    return(results_df)
  }

  # Apply the filtering and normality test for each dataset
  cat("\033[33mApply the filtering and normality test for each dataset\033[0m\n")
  imm_qti_filtered <- filter_valid_columns(imm_qti, design_formula_sym, "quanTIseq")
  imm_epic_filtered <- filter_valid_columns(imm_epic, design_formula_sym, "EPIC")
  imm_xcell_filtered <- filter_valid_columns(imm_xcell, design_formula_sym, "xcell")

  # Perform the normality test
  if (ncol(imm_qti_filtered) > 1) {
    cat("\033[32mPerforming Shapiro-Wilk test for quanTIseq\033[0m\n")
    resultados_norm_imm_qti <- prueba_norm(imm_qti_filtered, design_formula_sym)
    print(resultados_norm_imm_qti)
    # Perform parametric tests if normality is satisfied
    cat("\033[32mPerforming parametric tests for quanTIseq\033[0m\n")
    parametric_results_imm_qti <- perform_parametric_tests(imm_qti_filtered, design_formula_sym)
    print(parametric_results_imm_qti)
  } else {
    cat("\033[31mCan't perform Shapiro-Wilk test for quanTIseq. No valid columns available.\033[0m\n")
  }

  if (ncol(imm_epic_filtered) > 1) {
    cat("\033[32mPerforming Shapiro-Wilk test for EPIC\033[0m\n")
    resultados_norm_epic <- prueba_norm(imm_epic_filtered, design_formula_sym)
    print(resultados_norm_epic)
    # Perform parametric tests if normality is satisfied
    cat("\033[32mPerforming parametric tests for EPIC\033[0m\n")
    parametric_results_epic <- perform_parametric_tests(imm_epic_filtered, design_formula_sym)
    print(parametric_results_epic)
  } else {
    cat("\033[31mCan't perform Shapiro-Wilk test for EPIC. No valid columns available.\033[0m\n")
  }

  if (ncol(imm_xcell_filtered) > 1) {
    cat("\033[32mPerforming Shapiro-Wilk test for xcell\033[0m\n")
    resultados_norm_xcell <- prueba_norm(imm_xcell_filtered, design_formula_sym)
    print(resultados_norm_xcell)
    # Perform parametric tests if normality is satisfied
    cat("\033[32mPerforming parametric tests for xcell\033[0m\n")
    parametric_results_xcell <- perform_parametric_tests(imm_xcell_filtered, design_formula_sym)
    print(parametric_results_xcell)
  } else {
    cat("\033[31mCan't perform Shapiro-Wilk test for xcell. No valid columns available.\033[0m\n")
  }

  ############# Heatmaps
  # Transponer y estandarizar por filas
  h_imm_epic <- as.data.frame(t(imm_epic))
  ## Delete the row uncharacterized cell
  h_imm_epic <- head(h_imm_epic, -1)
  rn_himmepic <- rownames(h_imm_epic)
  # "Cancer associated fibroblast" will become CAFs
  rn_himmepic[2] <- "CAFs"
  cl_himmepic <- colnames(h_imm_epic)
  h_imm_epic <- apply(h_imm_epic, 2, as.numeric)
  h_imm_epic <- t(apply(h_imm_epic,1,scale))
  rownames(h_imm_epic) <- rn_himmepic
  colnames(h_imm_epic) <- cl_himmepic
  h_imm_epic <- as.data.frame(h_imm_epic)


  # Transponer y estandarizar por filas
  h_imm_qti <- as.data.frame(t(imm_qti))
  h_imm_qti <- head(h_imm_qti, -1)
  rn_himmqti <- rownames(h_imm_qti)
  rn_himmqti[7] <- "T cell CD4+"
  rn_himmqti[9] <- "T cell regulatory"
  cl_himmqti <- colnames(h_imm_qti)
  h_imm_qti <- apply(h_imm_qti, 2, as.numeric)
  h_imm_qti <- t(apply(h_imm_qti,1,scale))
  rownames(h_imm_qti) <- rn_himmqti
  colnames(h_imm_qti) <- cl_himmqti
  h_imm_qti <- as.data.frame(h_imm_qti)

  # Crear dataframe para heatmap
  h_imm_xcell <- as.data.frame(imm_xcell)
  # Transponer y estandarizar por filas
  h_imm_xcell <- as.data.frame(t(h_imm_xcell))
  h_imm_xcell <- head(h_imm_xcell, -4)

  # Poblaciones celulares no interesantes
  xcell_row_delete <- c('Common lymphoid progenitor', 'Common myeloid progenitor',
                        'Granulocyte-monocyte progenitor', 'Hematopoietic stem cell')

  h_imm_xcell <- dplyr::filter(h_imm_xcell, !rownames(h_imm_xcell) %in% xcell_row_delete)


  rn_himmxcell <- rownames(h_imm_xcell)
  cl_himmxcell <- colnames(h_imm_xcell)
  h_imm_xcell <- apply(h_imm_xcell, 2, as.numeric)
  h_imm_xcell <- t(apply(h_imm_xcell,1,scale))
  rownames(h_imm_xcell) <- rn_himmxcell
  colnames(h_imm_xcell) <- cl_himmxcell
  h_imm_xcell <- as.data.frame(h_imm_xcell)

  # Factorización de variables y renombrado en AnnotData
  AnnotData[[design_formula]] <- as.factor(AnnotData[[design_formula]])

  ################################### HEATMAP
  cat("\033[33mGENERATING HEATMAPS\033[0m\n")
  remove_nan_rows <- function(df) {
    # Remove rows with any NaN values
    df_cleaned <- df[!apply(df, 1, function(row) any(is.nan(row))), ]
    return(df_cleaned)
  }


  # Generar el heatmap para h_imm_epic
  combined_data <- remove_nan_rows(h_imm_epic)

  if (ncol(combined_data) < 3) {
    cat("Cannot generate the heatmap because the dataset has less than 3 columns after removing rows with NaN values.\n")
  } else {
    if (any(is.na(combined_data))) {
      cat("Cannot generate the heatmap because there are NaN values in the data.\n")
    } else {
      col_ann_data <- AnnotData[colnames(combined_data), , drop = FALSE]

      if (ncol(combined_data) == nrow(col_ann_data)) {
        col_ann_data <- as.data.frame(col_ann_data)
      } else {
        stop("Dimensions of col_ann_data and combined_data do not match")
      }

      annotation_col <- as.data.frame(col_ann_data[, design_formula, drop = FALSE])
      rownames(annotation_col) <- colnames(combined_data)

      HeatmapEPIC <- pheatmap::pheatmap(
        as.matrix(combined_data),
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        annotation_col = annotation_col,
        fontsize = 9,
        color = colorRampPalette(c("#4793AF", "white", "#013649"))(50),
        border_color = "grey60",
        main = "Heatmap EPIC",
        legend = TRUE,
        angle_col = 45,
        silent = TRUE
      )
    }
  }

  # Generar el heatmap para h_imm_qti
  combined_data <- remove_nan_rows(h_imm_qti)

  if (ncol(combined_data) < 3) {
    cat("Cannot generate the heatmap because the dataset has less than 3 columns after removing rows with NaN values.\n")
  } else {
    if (any(is.na(combined_data))) {
      cat("Cannot generate the heatmap because there are NaN values in the data.\n")
    } else {
      col_ann_data <- AnnotData[colnames(combined_data), , drop = FALSE]

      if (ncol(combined_data) == nrow(col_ann_data)) {
        col_ann_data <- as.data.frame(col_ann_data)
      } else {
        stop("Dimensions of col_ann_data and combined_data do not match")
      }

      annotation_col <- as.data.frame(col_ann_data[, design_formula, drop = FALSE])
      rownames(annotation_col) <- colnames(combined_data)

      Heatmap_qti <- pheatmap::pheatmap(
        as.matrix(combined_data),
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        annotation_col = annotation_col,
        fontsize = 9,
        color = colorRampPalette(c("#4793AF", "white", "#013649"))(50),
        border_color = "grey60",
        main = "Heatmap QTI",
        legend = TRUE,
        angle_col = 45,
        silent = TRUE
      )
    }
  }

  # Generar el heatmap para h_imm_xcell
  combined_data <- remove_nan_rows(h_imm_xcell)

  if (ncol(combined_data) < 3) {
    cat("Cannot generate the heatmap because the dataset has less than 3 columns after removing rows with NaN values.\n")
  } else {
    if (any(is.na(combined_data))) {
      cat("Cannot generate the heatmap because there are NaN values in the data.\n")
    } else {
      col_ann_data <- AnnotData[colnames(combined_data), , drop = FALSE]

      if (ncol(combined_data) == nrow(col_ann_data)) {
        col_ann_data <- as.data.frame(col_ann_data)
      } else {
        stop("Dimensions of col_ann_data and combined_data do not match")
      }

      annotation_col <- as.data.frame(col_ann_data[, design_formula, drop = FALSE])
      rownames(annotation_col) <- colnames(combined_data)

      Heatmap_xcell <- pheatmap::pheatmap(
        as.matrix(combined_data),
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        annotation_col = annotation_col,
        fontsize = 9,
        color = colorRampPalette(c("#4793AF", "white", "#013649"))(50),
        border_color = "grey60",
        main = "Heatmap XCELL",
        legend = TRUE,
        angle_col = 45,
        silent = TRUE
      )
    }
  }

  cat("\033[32mHeatmap of qti, EPIC and xcell will be stored on plots_TME_heatmap.pdf\033[0m\n")

  if (exists("Heatmap_qti")) {
    pdf("plots_TME_Heatmap_qti.pdf", width = 11, height = 14)
    print(Heatmap_qti)
    dev.off()
  } else {
    cat("Heatmap_qti object does not exist.\n")
  }

  if (exists("HeatmapEPIC")) {
    pdf("plots_TME_HeatmapEPIC.pdf", width = 11, height = 14)
    print(HeatmapEPIC)
    dev.off()
  } else {
    cat("HeatmapEPIC object does not exist.\n")
  }
  if (exists("Heatmap_xcell")) {
    pdf("plots_TME_Heatmap_xcell.pdf", width = 11, height = 14)
    print(Heatmap_xcell)
    dev.off()
  } else {
    cat("Heatmap_xcell object does not exist.\n")
  }

  data_table_list <- list(EPIC = imm_epic, QTI = imm_qti, XCELL = imm_xcell)
  return(data_table_list)
}
