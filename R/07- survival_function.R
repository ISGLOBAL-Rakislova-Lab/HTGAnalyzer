#' Perform Survival Analysis
#'
#' This function performs survival analysis based on the provided data. The function uses parameters in a hierarchical manner:
#' if `res` is provided, it uses the top genes from `res`; if `genes_to_use` is provided, it uses those specific genes;
#' if `TME` is provided, it uses the genes from `TME`. If multiple parameters are provided, `res` takes precedence over
#' `genes_to_use`, and `genes_to_use` takes precedence over `TME`.
#'
#' @param variable_01 Character. The name of the survival event variable (e.g., "Recurrence_01").
#' @param time Character. The name of the time variable (e.g., "Time_to_death_surv").
#' @param col_data Data frame. Annotation data containing the survival event and time variables.
#' @param count_data Matrix. Raw count data for the genes.
#' @param DEA DESeqDataSet object or NULL. Pre-existing DESeqDataSet object. If NULL, a new DESeqDataSet will be created from `count_data` and `col_data`.
#' @param res Data frame or NULL. Result data with `padj` values to identify top genes. Used if `genes_to_use` is NULL.
#' @param genes_to_use Character vector or NULL. Specific genes to use for survival analysis. Takes precedence over `res`.
#' @param TME Data frame or NULL. Data with gene expressions used if neither `res` nor `genes_to_use` are provided.
#' @param outliers Character vector. Outlier identifiers to be removed from the analysis.
#' @param pattern Character. Pattern to match and remove certain rows from the count data (e.g., "^NC-|^POS-|^GDNA-|^ERCC-").
#' @param remove_outliers Logical. Indicates whether outliers should be removed. Defaults to TRUE.
#'
#' @return Generates Kaplan-Meier survival plots and saves them as PDF files. The plots are saved in the current working directory.
#' @examples
#' HTG_survival(variable_01 = "Recurrence_01",
#'              time = "Time_to_death_surv",
#'              col_data = AnnotData,
#'              count_data = counts_data,
#'              res = res,
#'              genes_to_use = NULL,
#'              outliers = outliers,
#'              pattern = "^NC-|^POS-|^GDNA-|^ERCC-",
#'              remove_outliers = TRUE)
#'
#' HTG_survival(variable_01 = "Recurrence_01",
#'              time = "Time_to_death_surv",
#'              col_data = AnnotData,
#'              count_data = counts_data,
#'              res = res,
#'              genes_to_use = c("LCP1", "OMA1"),
#'              outliers = outliers,
#'              pattern = "^NC-|^POS-|^GDNA-|^ERCC-",
#'              remove_outliers = TRUE)
#'
#' HTG_survival(variable_01 = "Recurrence_01",
#'              time = "Time_to_death_surv",
#'              col_data = AnnotData,
#'              count_data = counts_data,
#'              res = NULL,
#'              genes_to_use = NULL,
#'              TME = TME_data$EPIC,
#'              outliers = outliers,
#'              pattern = "^NC-|^POS-|^GDNA-|^ERCC-",
#'              remove_outliers = TRUE)
#'
#' @export
HTG_survival <- function(variable_01, time, col_data, count_data, DEA = NULL, res = NULL, genes_to_use = NULL,
                         outliers, pattern = NULL, remove_outliers = TRUE, TME = NULL) {
  library(dplyr)
  library(survival)
  library(maxstat)
  library(survminer)
  library(ggplot2)
  library(gridExtra)


  if (is.null(variable_01) || is.null(time)) {
    stop("Variables for survival analysis are required.")
  }

  if (!is.null(pattern)) {
    # Remove outliers based on pattern
    filtered <- subset(count_data, !grepl(pattern, rownames(count_data)))
  } else {
    filtered <- count_data
  }

    if (remove_outliers) {
    counts_filtered <- filtered[, !colnames(filtered) %in% outliers]
    AnnotData <- col_data[!col_data[["id"]] %in% outliers, ]
  } else {
    counts_filtered <- count_data
    AnnotData <- col_data
  }

  cat("\033[32mNormalizing data...\033[0m\n")
  design_formul <- as.formula("~ 1")
  colnames(AnnotData) <- gsub(" ", "_", colnames(AnnotData))
  col_data <- AnnotData[order(AnnotData$id), ]
  col_data<- as.data.frame(col_data)
  counts_filtered <- counts_filtered[, order(colnames(counts_filtered))]
  if (!identical(colnames(counts_filtered), col_data$id)) {
    stop("Column names of counts_filtered and IDs in col_data do not match.")
  }

  # Clean column names to avoid issues with special characters
  clean_column_names <- function(names) {
    gsub("[^[:alnum:]_]", "_", names)
  }

  # Create DESeqDataSet object
  if (is.null(DEA)) {
    # Create DESeqDataSet object
    rownames(col_data) <- col_data$id
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts_filtered, colData = col_data, design = ~ 1)
    dds <- DESeq2::estimateSizeFactors(dds)
  }
  normalized_counts <- DESeq2::counts(dds, normalized = TRUE)
  df_t <- t(normalized_counts)

  # Checking for duplicate columns
  cat("\033[32mChecking for duplicate columns\033[0m\n")
  duplicated_columns <- colnames(df_t)[duplicated(colnames(df_t))]
  if (length(duplicated_columns) > 0) {
    print(paste("Duplicate columns found:", paste(duplicated_columns, collapse = ", ")))
  } else {
    print("No duplicate columns.")
  }

  rownames(col_data) <- col_data$id
  ids_data <- rownames(df_t)


  # Filter data
  cat("\033[32mFiltering data.\033[0m\n")
  subset_data <- col_data %>% dplyr::filter(id %in% ids_data)

  df_ta <- as.data.frame(df_t)
  df_ta$id <- rownames(df_ta)



  # Selecting genes
  if (!is.null(res)) {
    cat("\033[32mSelecting TOP 10 genes with the lowest padj\033[0m\n")
    top_genes <- rownames(head(res[order(res$padj), ], 10))
    selected_df_t <- df_t[, top_genes, drop = FALSE]
    selected_df_t <- as.data.frame(selected_df_t)
    selected_df_t$id <- rownames(selected_df_t)
    col_data$id <- rownames(col_data)
    merged_data <- merge(col_data, selected_df_t, by = "id")
    merged_data[[time]] <- as.numeric(merged_data[[time]])

    cat("\033[32mStarting survival analysis.\033[0m\n")

    pdf("survival_analysis_plots.pdf")

    # Replace special characters in top_genes
    top_genes_clean <- clean_column_names(top_genes)

    # Perform survival analysis for each gene
    for (i in top_genes_clean) {
      if (!is.numeric(merged_data[[i]])) {
        next
      }

      cat("\n")
      cat("\033[32mPerforming analysis for column:\033[0m ", i, "\n")

      # Perform MAXSTAT test
      merged_data$time <- merged_data[[time]]
      merged_data$variable_01 <- merged_data[[variable_01]]
      gene_column <- get(i, merged_data)  # Access the column dynamically
      MAXSTAT <- maxstat.test(Surv(time, variable_01) ~ gene_column, data = merged_data,
                              smethod = "LogRank", pmethod = "Lau92", iscores = TRUE, minprop = 0.45, maxprop = 0.55)
      cut.off <- MAXSTAT$estimate
      cat("\033[32mCUT OFF\033[0m\n")
      print(cut.off)

      # Create a new variable based on the cutoff
      new_column_name <- paste0(i, "_mRNA_expression")
      merged_data[[new_column_name]] <- ifelse(gene_column > cut.off, "High", "Low")
      merged_data[[new_column_name]] <- factor(merged_data[[new_column_name]])

      # Fit survival model
      cat("\033[32mFitting survival model\033[0m\n")
      surv_object <- Surv(merged_data$time, merged_data$variable_01)
      surv_formula <- as.formula(paste("surv_object ~", new_column_name))

      fit1 <- survfit(surv_formula, data = merged_data)

      # Summary of the fit
      cat("\033[32mSummary of the fit\033[0m\n")
      print(summary(fit1))

      # Log-rank test and p-value
      cat("\033[32mPerforming log-rank test and obtaining p-value\033[0m\n")
      surv_diff <- survdiff(surv_formula, data = merged_data)
      p_value <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)

      print(surv_diff)
      cat("\033[32mP-value\033[0m\n")
      print(p_value)

      # Generate Kaplan-Meier plot
      cat("\033[32mGenerating Kaplan-Meier plot\033[0m\n")
      palette <- c("#9A3449", "#D4A8B1")
      plot(fit1, lty = 1, col = palette, lwd = 4, main = paste("Survival analysis for", i, "\n", "p-value =", format(p_value, digits = 3)))

      # Add a legend
      legend("topright",
             legend = c("High", "Low"),
             lty = 1,
             col = palette,
             lwd = 4)
      cat("\033[32mPlots saved in survival_analysis_plots.pdf\033[0m\n")
    }
    dev.off()
    while (!is.null(dev.list())) dev.off()

  } else if (!is.null(genes_to_use)) {
    cat("\033[32mUsing provided genes\033[0m\n")
    top_genes <- genes_to_use
    selected_df_t <- df_t[, top_genes, drop = FALSE]
    selected_df_t <- as.data.frame(selected_df_t)
    selected_df_t$id <- rownames(selected_df_t)
    col_data$id <- rownames(col_data)
    merged_data <- merge(col_data, selected_df_t, by = "id")
    merged_data[[time]] <- as.numeric(merged_data[[time]])


    cat("\033[32mStarting survival analysis.\033[0m\n")

    pdf("survival_analysis_plots_SELECTED_GENES.pdf")

    # Replace special characters in top_genes
    top_genes_clean <- clean_column_names(top_genes)

    # Perform survival analysis for each gene
    for (i in top_genes_clean) {
      if (!is.numeric(merged_data[[i]])) {
        next
      }

      cat("\n")
      cat("\033[32mPerforming analysis for column:\033[0m ", i, "\n")

      # Perform MAXSTAT test
      merged_data$time <- merged_data[[time]]
      merged_data$variable_01 <- merged_data[[variable_01]]
      gene_column <- get(i, merged_data)  # Access the column dynamically
      MAXSTAT <- maxstat.test(Surv(time, variable_01) ~ gene_column, data = merged_data,
                              smethod = "LogRank", pmethod = "Lau92", iscores = TRUE, minprop = 0.45, maxprop = 0.55)
      cut.off <- MAXSTAT$estimate
      cat("\033[32mCUT OFF\033[0m\n")
      print(cut.off)

      # Create a new variable based on the cutoff
      new_column_name <- paste0(i, "_mRNA_expression")
      merged_data[[new_column_name]] <- ifelse(gene_column > cut.off, "High", "Low")
      merged_data[[new_column_name]] <- factor(merged_data[[new_column_name]])

      # Fit survival model
      cat("\033[32mFitting survival model\033[0m\n")
      surv_object <- Surv(merged_data$time, merged_data$variable_01)
      surv_formula <- as.formula(paste("surv_object ~", new_column_name))

      fit1 <- survfit(surv_formula, data = merged_data)

      # Summary of the fit
      cat("\033[32mSummary of the fit\033[0m\n")
      print(summary(fit1))

      # Log-rank test and p-value
      cat("\033[32mPerforming log-rank test and obtaining p-value\033[0m\n")
      surv_diff <- survdiff(surv_formula, data = merged_data)
      p_value <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)

      print(surv_diff)
      cat("\033[32mP-value\033[0m\n")
      print(p_value)

      # Generate Kaplan-Meier plot
      cat("\033[32mGenerating Kaplan-Meier plot\033[0m\n")
      palette <- c("#9A3449", "#D4A8B1")
      plot(fit1, lty = 1, col = palette, lwd = 4, main = paste("Survival analysis for", i, "\n", "p-value =", format(p_value, digits = 3)))

      # Add a legend
      legend("topright",
             legend = c("High", "Low"),
             lty = 1,
             col = palette,
             lwd = 4)
      cat("\033[32mPlots saved in survival_analysis_plots.pdf\033[0m\n")
    }
    dev.off()
    while (!is.null(dev.list())) dev.off()

  } else if (!is.null(TME)) {
    cat("\033[32mUsing rownames from TME\033[0m\n")
    TME <- TME[, -ncol(TME)]
    colnames(TME) <- gsub(" ", "_", colnames(TME))
    col_data$id <- rownames(col_data)
    TME$id <- rownames(TME)
    merged_data <- merge(col_data, TME, by = "id")
    merged_data[[time]] <- as.numeric(merged_data[[time]])

    cat("\033[32mStarting survival analysis.\033[0m\n")
    colnames(merged_data) <- gsub(" ", "_", colnames(merged_data))
    colnames(merged_data) <- gsub("\\+", "_", colnames(merged_data))


    pdf("survival_analysis_plots_TME.pdf")

    for (i in colnames(TME)) {
      if (!is.numeric(merged_data[[i]])) {
        next
      }

      cat("\n")
      cat("\033[32mPerforming analysis for column:\033[0m ", i, "\n")

      # Perform MAXSTAT test
      merged_data$time<- merged_data[[time]]
      merged_data$variable_01<- merged_data[[variable_01]]
      MAXSTAT <- maxstat.test(Surv(time, variable_01) ~ merged_data[[i]], data = merged_data,
                              smethod = "LogRank", pmethod = "Lau92", iscores = TRUE, minprop = 0.45, maxprop = 0.55)
      cut.off <- MAXSTAT$estimate
      cat("\033[32mCUT OFF\033[0m\n")
      print(cut.off)

      # Create a new variable based on the cutoff

      merged_data[[paste0(i)]] <- ifelse(merged_data[[i]] > cut.off, "High", "Low")
      merged_data[[paste0(i)]] <- factor(merged_data[[paste0(i)]])


      # Fit survival model
      cat("\033[32mFitting survival model\033[0m\n")
      column_name <- paste0(i)
      surv_object <- Surv( merged_data$time , merged_data$variable_01)
      surv_formula <- as.formula(paste("surv_object ~", column_name))

      fit1 <- survfit(surv_formula, data = merged_data)
      # Summary of the fit

      cat("\033[32mSummary of the fit\033[0m\n")
      print(summary(fit1))

      # Log-rank test and p-value
      cat("\033[32mPerforming log-rank test and obtaining p-value\033[0m\n")
      surv_diff <- survdiff(surv_formula, data = merged_data)
      p_value <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)

      print(surv_diff)
      cat("\033[32mP-value\033[0m\n")
      print(p_value)

      # Generate Kaplan-Meier plot
      cat("\033[32mGenerating Kaplan-Meier plot\033[0m\n")
      palette <- c("#9A3449", "#D4A8B1")
      plot(fit1, lty = 1, col = palette, lwd = 4, main = paste("Survival analysis for", i, "\n", "p-value =", format(p_value, digits = 3)))

      # Añadir una leyenda
      legend("topright",
             legend = c("High", "Low"),
             lty = 1,
             col = palette,
             lwd = 4)
      cat("\033[32mPlots saved in survival_analysis_plots.pdf\033[0m\n")
    }
    dev.off()
    while (!is.null(dev.list())) dev.off()


  } else {
    stop("Either 'res', 'genes_to_use', or 'TME' must be provided.")
  }
}
