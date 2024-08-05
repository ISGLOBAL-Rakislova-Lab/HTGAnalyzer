#' HTG_TME
#'
#' @description This function carries out a comprehensive convolution analysis involving TPM normalization and deconvolution using multiple methods (EPIC, quanTIseq, and xCell). It produces several output files containing the results of the normalization and deconvolution processes.

#' @param outliers Vector of outlier sample IDs to be removed from analysis.
#' @param pattern Regular expression pattern to identify outlier sample IDs.
#' @param counts_filtered Count data matrix of gene expression.
#' @param AnnotData DataFrame containing sample metadata.
#' @param design_formula Formula specifying the design of the analysis.
#' @param correlation Logical; whether to perform correlation analysis.
#' @param dds Preprocessed DESeqDataSet for differential expression analysis.
#' @param generate_volcano Logical; whether to generate volcano plots.
#' @param remove_outliers Logical; whether to remove outlier samples.
#'
#' @return The function outputs multiple CSV files:
#' - `tpm_counts.csv`: TPM normalized counts.
#' - `imm_epic.csv`, `imm_qti.csv`, `imm_xcell.csv`: Results of deconvolution methods EPIC, quanTIseq, and xCell respectively.
#'
#' @export
#' @examples
#' HTG_TME(outliers, pattern= "^NC-|^POS-|^GDNA-|^ERCC-",
#'                  counts_filtered, AnnotData, design_formula= "Ciclina2", remove_outliers = TRUE, dds = NULL)
#'
#' @name HTG_TME
#'
#'
#'
HTG_TME <- function(outliers, pattern = NULL, counts_filtered, AnnotData, design_formula = NULL ,
                            correlation = FALSE, dds = NULL, generate_volcano = TRUE, remove_outliers = TRUE) {
  suppressMessages({
    library(DESeq2)
    library(ggplot2)
    library(ggrepel)
    library(PoiClaClu)
    library(EnhancedVolcano)
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(msigdbr)
    library(fgsea)
    library(DOSE)
    library(enrichplot)
    library(pheatmap)
    library(RColorBrewer)
    library(ggupset)
    library(grid)
    library(survival)
    library(survminer)
    library(dplyr)
    library(IOBR)
    library(immunedeconv)
    library(tidyverse)
    library(ggpubr)
    library(tidyr)
    library(maxstat)
  })

  if (remove_outliers) {
    if (!is.null(pattern)) {
      # Remove outliers based on pattern
      filtered <- subset(counts_filtered, !grepl(pattern, rownames(counts_filtered)))
    } else {
      # If no pattern is provided, do not filter based on pattern
      filtered <- counts_filtered
    }
    # Remove columns corresponding to outliers
    counts_filtered <- filtered[, !colnames(filtered) %in% outliers]
    AnnotData <- AnnotData[!AnnotData[["id"]] %in% outliers, ]
  } else {
    counts_filtered <- counts_filtered
    AnnotData <- AnnotData
  }

  cat("\033[32mCOnvolution analysis skipped.\033[0m\n")
  print(paste0("Número inicial de genes: ", dim(counts_filtered)[1]))
  ## Normalización TPM
  cat("\033[32mTPM normalization performed and stored on tpm_counts.csv\033[0m\n")
  suppressWarnings({
    tpm_counts <- count2tpm(counts_filtered,
                            idType = "Symbol",
                            org = "hsa",
                            source = "biomart")
    write.csv(tpm_counts, "tpm_counts.csv")
  })

  # Se almacenan los genes omitidos en la normalización TPM
  genes_omitidos <- setdiff(rownames(counts_filtered), rownames(tpm_counts))
  print(paste0("Number of genes omitted during TPM normalization due to their length not being available in Biomart:   ", dim(counts_filtered)[1]-dim(tpm_counts)[1]))
  tpm_counts <- as.data.frame(tpm_counts)

  if (!is.null(dds)) {
    cat("\033[32mWe are going to use information from DEA\033[0m\n")
    dds <- estimateSizeFactors(dds)
    normalized_counts <- counts(dds, normalized = TRUE)
  } else {
    cat("Performing normalization.\n")
    design_formul <- as.formula("~ 1")
    colnames(AnnotData) <- gsub(" ", "_", colnames(AnnotData))
    AnnotData <- AnnotData[order(AnnotData$id), ]
    counts_filtered <- counts_filtered[, order(colnames(counts_filtered))]
    if (!identical(colnames(counts_filtered), AnnotData$id)) {
      stop("Column names of counts_filtered and IDs in AnnotData do not match.")
    }
  }

  ## Deconvolución
  imm_epic <- deconvolute(tpm_counts, method = "epic")
  imm_qti <- deconvolute(tpm_counts, method = "quantiseq")
  imm_xcell <- deconvolute(tpm_counts, method = "xcell")
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
  dim(AnnotData)
  dim(imm_epic)

  col_data_ids <- AnnotData$id
  imm_epic_ids <- as.character(rownames(imm_epic))
  common_ids <- intersect(col_data_ids, imm_epic_ids)
  # Filtrar los datos originales según los IDs comunes
  AnnotData <- AnnotData[col_data_ids %in% common_ids, ]
  dim(AnnotData)
  AnnotData<- as.data.frame(AnnotData)
  rownames(AnnotData)<- AnnotData$id
  AnnotData <- AnnotData[order(rownames(AnnotData)), ]
  imm_epic <- imm_epic[order(rownames(imm_epic)), ]
  imm_qti <- imm_qti[order(rownames(imm_qti)), ]
  imm_xcell <- imm_xcell[order(rownames(imm_xcell)), ]

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

  ############# imm_epic
  design_formula_sym <- sym(design_formula)
  # Vector with column names except the last one
  column_names <- names(imm_epic)[-ncol(imm_epic)]
  # Loop to generate plots for each column
  for (col_name in column_names) {
    # Calculate means by design_formula
    means_df <- imm_epic %>%
      group_by(!!design_formula_sym) %>%
      summarize(mean_value = mean(.data[[col_name]] * 100, na.rm = TRUE), .groups = "drop")
    # Generate plot
    print(
      ggplot(imm_epic, aes(x = !!design_formula_sym , y = .data[[col_name]] * 100,
                           fill = !!design_formula_sym, color = !!design_formula_sym)) +
        geom_jitter(alpha = 1, width = 0.3, height = 0) +
        geom_boxplot(fill = "white", alpha = 0.5, outlier.alpha = 1) +
        geom_point(data = means_df, aes(x = !!design_formula_sym, y = mean_value, fill = !!design_formula_sym),
                   shape = 22, color = "black", size = 3, stroke = 1.5,
                   show.legend = F) +
        labs(x = NULL, y = "Abundance (%)", title = col_name, subtitle = "imm_epic") +
        scale_y_continuous(labels = scales::percent_format(scale = 1),
                           limits = c(0, NA)) +
        theme(axis.text.x = element_blank())
    )
  }
  ############# imm_qti
  # Vector with column names except the last one
  column_names <- names(imm_qti)[-ncol(imm_qti)]
  # Loop to generate plots for each column
  pdf("plots_imm_qti.pdf")
  for (col_name in column_names) {
    # Calculate means by design_formula
    means_df <- imm_qti %>%
      group_by(!!design_formula_sym) %>%
      summarize(mean_value = mean(.data[[col_name]] * 100, na.rm = TRUE))
    # Generate plot
    print(
      ggplot(imm_qti, aes(x = !!design_formula_sym , y = .data[[col_name]] * 100,
                          fill = !!design_formula_sym, color = !!design_formula_sym)) +
        geom_jitter(alpha = 1, width = 0.3, height = 0) +
        geom_boxplot(fill = "white", alpha = 0.5, outlier.alpha = 1) +
        geom_point(data = means_df, aes(x = !!design_formula_sym, y = mean_value, fill = !!design_formula_sym),
                   shape = 22, color = "black", size = 3, stroke = 1.5,
                   show.legend = F) +
        labs(x = NULL, y = "Abundance (%)",  title = col_name, subtitle = "imm_qti") +
        scale_y_continuous(labels = scales::percent_format(scale = 1),
                           limits = c(0, NA)) +
        theme(axis.text.x = element_blank())
    )
  }
  dev.off()

  ############# imm_xcell
  # Vector con los nombres de las columnas excepto la última
  column_names <- names(imm_xcell)[-ncol(imm_xcell)]
  # Bucle for para generar gráficos para cada columna
  for (col_name in column_names) {
    # Calcular las medias por design_formula
    means_df <- imm_xcell %>%
      group_by(!!design_formula_sym) %>%
      summarize(mean_value = mean(.data[[col_name]] * 100, na.rm = TRUE))

    print(
      ggplot(imm_xcell, aes(x = !!design_formula_sym , y = .data[[col_name]],
                            fill = !!design_formula_sym, color = !!design_formula_sym)) +
        geom_jitter(alpha = 1, width = 0.3, height = 0) +
        geom_boxplot(fill = "white", alpha = 0.5, outlier.alpha = 1) +
        geom_point(data = means_df, aes(x = !!design_formula_sym, y = mean_value, fill = !!design_formula_sym),
                   shape = 22, color = "black", size = 3, stroke = 1.5,
                   show.legend = F) +
        labs(x = NULL, y = "Abundance (arbitrary value)", title = col_name, subtitle = "imm_xcell") +
        scale_y_continuous(expand = expansion(add = c(0, 0.1)))  +
        theme(axis.text.x = element_blank())
    )
  }


  # Composición celular del TME i por grupo

  # Function to generate the bar plot for each sample
  plot_bar <- function(df, paleta, titulo, legend.position) {
    df <- df %>%
      rownames_to_column(var = "Sample") %>%
      pivot_longer(cols = colnames(df)[2:(ncol(df) - 1)],
                   names_to = "Cell_Type", values_to = "Value") %>%
      mutate(Sample = factor(Sample, levels = rev(unique(Sample))),
             Cell_Type = factor(Cell_Type, levels = rev(unique(Cell_Type))))

    p <- ggplot(df, aes(x = Sample, y = Value, fill = Cell_Type)) +
      geom_bar(stat = "identity") +
      labs(title = titulo,
           x = "Samples",
           y = "Cell Fraction (%)") +
      coord_flip() +
      guides(fill = guide_legend(reverse = TRUE)) +
      scale_fill_manual(values = paleta) +
      theme_minimal() +
      theme(legend.position = legend.position,
            axis.text.y = element_text(size = 4)) +
      scale_y_continuous(labels = scales::percent)

    return(p)
  }

  # Function to generate the average bar plot by group (e.g., HPV status)
  plot_bar_group <- function(df, paleta, titulo, design_formula, legend_position = "right") {
    design_formula_sym <- sym(design_formula)
    niveles_tipo_cel <- colnames(df)[1:(ncol(df) - 1)]

    promedios <- df %>%
      rownames_to_column(var = "Sample") %>%
      pivot_longer(cols = niveles_tipo_cel,
                   names_to = "Cell_Type", values_to = "Value") %>%
      group_by(!!design_formula_sym, Cell_Type) %>%
      summarise(Average = mean(Value, na.rm = TRUE)) %>%
      ungroup() %>%
      mutate(!!design_formula_sym := factor(!!design_formula_sym, levels = rev(unique(!!design_formula_sym))),
             Cell_Type = factor(Cell_Type, levels = rev(niveles_tipo_cel)))

    p <- ggplot(promedios, aes(x = !!design_formula_sym, y = Average, fill = Cell_Type)) +
      geom_bar(stat = "identity") +
      labs(title = titulo,
           x = design_formula,
           y = "Average Cell Fraction (%)") +
      coord_flip() +
      guides(fill = guide_legend(reverse = TRUE)) +
      scale_fill_manual(values = paleta) +
      theme_minimal() +
      theme(legend.position = legend_position,
            axis.text.y = element_text(size = 8)) +
      scale_y_continuous(labels = scales::percent)

    return(p)
  }

  # Function to combine both plots into one
  plot_combined <- function(df, paleta, titulo_individual, titulo_grupo, design_formula, legend_position = "right") {
    p1 <- plot_bar(df, paleta, titulo_individual, legend_position)
    p2 <- plot_bar_group(df, paleta, titulo_grupo, design_formula, legend_position)


    combined_plot <- ggarrange(p1, p2, ncol = 1, nrow = 2, heights = c(1, 1))
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
  print(combined_plot_EPIC)

  combined_plot_quanTIseq <- plot_combined(imm_qti, paleta_qti, "quanTIseq Individual", "quanTIseq Average", design_formula, "right")
  print(combined_plot_quanTIseq)

  combined_plot_xCell <- plot_combined(imm_xcell, paleta_extendida, "xCell Individual", "xCell Average", design_formula, "right")
  print(combined_plot_xCell)

  pdf("plot_cell_fraction_Average_cell_fraction.pdf")
  print(combined_plot_EPIC)
  print(combined_plot_quanTIseq)
  print(combined_plot_xCell)
  dev.off()

  # # Contrastes de Hipótesis
  prueba_norm <- function(df, design_formula_sym) {
    # Transform the dataframe to long format
    df_largo <- df %>%
      pivot_longer(-design_formula_sym,
                   names_to = "Celular type", values_to = "Fraction")

    # Perform Shapiro-Wilk test for each cell type
    resultados_normalidad <- df_largo %>%
      group_by(`Celular type`) %>%
      summarise(shapiro_test = list(shapiro.test(Fraction))) %>%
      mutate(p.value = map_dbl(shapiro_test, "p.value"))

    return(resultados_normalidad)
  }
  check_column_equal <- function(column) {
    all_equal <- length(unique(column)) == 1
    return(all_equal)
  }

  # Verificación de igualdad para cada columna
  cat("\033[32mChecking if All Values in Each Column Are Equal\033[0m\n")
  equality_results_qti <- sapply(imm_qti, check_column_equal)
  unequal_columns_qti <-sum(equality_results_qti)

  equality_results_epic <- sapply(imm_epic, check_column_equal)
  unequal_columns_epic<- sum(equality_results_epic)

  equality_results_xcell <- sapply(imm_xcell, check_column_equal)
  print(equality_results_xcell)
  unequal_columns_xcell<-sum(equality_results_xcell)


  if (unequal_columns_qti == 0) {
    resultados_norm_imm_qti <- prueba_norm(imm_qti, design_formula_sym)
    cat("\033[32mTest for normality: Shapiro-Wilk test of quanTIseq\033[0m\n")
    print(resultados_norm_imm_qti)
  }else{
    cat("\033[32mCan't perform Shapiro-Wilk test for . Some columns have same value\033[0m\n")
    print(equality_results_qti)
  }


  if (unequal_columns_epic == 0) {
    resultados_norm_epic <- prueba_norm(imm_epic, design_formula_sym)
    cat("\033[32mTest for normality: Shapiro-Wilk test of EPIC\033[0m\n")
    print(resultados_norm_epic)
  }else{
    cat("\033[32mCan't perform Shapiro-Wilk test fro EPIC. Some columns have same value")
    print(equality_results_epic)
  }
  if (unequal_columns_xcell == 0) {
    resultados_norm_xcell <- prueba_norm(imm_xcell, design_formula_sym)
    cat("\033[32mTest for normality: Shapiro-Wilk test of xcell\033[0m\n")
    print(resultados_norm_xcell)
  }else{
    cat("\033[32mCan't perform Shapiro-Wilk test fro xcell. Some columns have same value\033[0m\n")
    print(equality_results_epic)
  }

  # Heatmaps
  # Function to transform dataframe to long format
  trans_formato_largo <- function(df, design_formula_sym) {
    # df_largo <- df %>%
    #   pivot_longer(-design_formula_sym, names_to = "Celular type", values_to = "Fraction")
    df_largo <- df %>%
      pivot_longer(-all_of(design_formula_sym),
                   names_to = "Celular type", values_to = "Fraction")
    return(df_largo)
  }
  # Example usage to transform dataframes to long format
  imm_epic<- as.data.frame(imm_epic)
  imm_qti<- as.data.frame(imm_qti)
  imm_qti<- as.data.frame(imm_qti)
  epic_largo <- trans_formato_largo(imm_epic, design_formula_sym)
  qti_largo <- trans_formato_largo(imm_qti, design_formula_sym)
  xcell_largo <- trans_formato_largo(imm_xcell, design_formula_sym)

  ## Comparación entre Grupos
  perform_parametric_test <- function(df, variable, groups) {
    resultados <- list()
    if (length(groups) != 2) {
      stop("groups should be a list with exactly two elements.")
    }
    # Iterate over column names except the last one (assuming the last one  or a similar grouping variable)
    for (col_name in colnames(df)[1:(ncol(df) - 1)]) {
      # Perform t-test for group 1
      if (!is.null(groups[[1]]) && length(groups[[1]]) > 1) {
        resultado <- t.test(df[[col_name]] ~ df[[variable]], subset = df[[variable]] %in% groups[[1]])
        resultados[[paste0(col_name, "_", groups[[1]][1], "_vs_", groups[[1]][2])]] <- c(col_name, paste0(groups[[1]][1], " vs ", groups[[1]][2]), resultado$p.value)
      } else {
        warning("Group 1 does not have enough samples.")
      }

      # Perform t-test for group 2
      if (!is.null(groups[[2]]) && length(groups[[2]]) > 1) {
        resultado <- t.test(df[[col_name]] ~ df[[variable]], subset = df[[variable]] %in% groups[[2]])
        resultados[[paste0(col_name, "_", groups[[2]][1], "_vs_", groups[[2]][2])]] <- c(col_name, paste0(groups[[2]][1], " vs ", groups[[2]][2]), resultado$p.value)
      } else {
        warning("Group 2 does not have enough samples.")
      }
    }
    resultados_df <- do.call(rbind, resultados)
    colnames(resultados_df) <- c("Tipo_Celular", "Comparacion", "P-valor")
    rownames(resultados_df) <- NULL
    resultados_df <- as.data.frame(resultados_df)
    print(resultados_df)
  }

  # Heatmaps
  # Crear dataframe para heatmap
  h_imm_epic <- imm_epic
  # Transponer y estandarizar por filas
  h_imm_epic <- as.data.frame(t(h_imm_epic))
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

  # Crear dataframe para heatmap
  h_imm_qti <- imm_qti
  # Transponer y estandarizar por filas
  h_imm_qti <- as.data.frame(t(h_imm_qti))
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
  # Renombrar nombres de poblaciones celulares antes de transponer
  h_imm_xcell <- h_imm_xcell %>%
    dplyr::rename('Myeloid DC activated' = 'Myeloid dendritic cell activated') %>%
    dplyr::rename('T cell CD4+' = 'T cell CD4+ (non-regulatory)') %>%
    dplyr::rename('Myeloid DC' = 'Myeloid dendritic cell') %>%
    dplyr::rename('CAFs' = 'Cancer associated fibroblast') %>%
    dplyr::rename('Plasmacytoid DC' = 'Plasmacytoid dendritic cell') %>%
    dplyr::rename('T cell regulatory' = 'T cell regulatory (Tregs)')

  # Transponer y estandarizar por filas
  h_imm_xcell <- as.data.frame(t(h_imm_xcell))
  h_imm_xcell <- head(h_imm_xcell, -4)

  # Poblaciones celulares no interesantes
  xcell_row_delete <- c('Common lymphoid progenitor', 'Common myeloid progenitor',
                        'Granulocyte-monocyte progenitor', 'Hematopoietic stem cell')

  h_imm_xcell <- h_imm_xcell %>%
    filter(!row.names(.) %in% xcell_row_delete)  # Eliminar filas

  rn_himmxcell <- rownames(h_imm_xcell)
  cl_himmxcell <- colnames(h_imm_xcell)
  h_imm_xcell <- apply(h_imm_xcell, 2, as.numeric)
  h_imm_xcell <- t(apply(h_imm_xcell,1,scale))
  rownames(h_imm_xcell) <- rn_himmxcell
  colnames(h_imm_xcell) <- cl_himmxcell
  h_imm_xcell <- as.data.frame(h_imm_xcell)

  # Factorización de variables y renombrado en AnnotData
  AnnotData[[design_formula]] <- as.factor(AnnotData[[design_formula]])

  ###################################
  #### EPIC
  combined_data <- h_imm_epic
  # Verificar si hay NaN en combined_data
  if (any(is.na(combined_data))) {
    cat("Cannot generate the EPIc heatmap because there are NaN values in the data.\n")
  } else {
    # Continuar con la generación del heatmap
    combined_data <- h_imm_epic
    col_ann_data <- AnnotData[colnames(combined_data), , drop = FALSE]

    # Asegúrate de que col_ann_data tenga las mismas dimensiones que combined_data
    if (ncol(combined_data) == nrow(col_ann_data)) {
      col_ann_data <- as.matrix(col_ann_data)  # Convertir a matriz
    } else {
      stop("Dimensions of col_ann_data and combined_data do not match")
    }

    # Anotaciones
    col_ann <- HeatmapAnnotation(
      df = col_ann_data[, design_formula, drop = FALSE],  # Selecciona solo la columna design_formula
      gp = gpar(col = "grey60", lwd = 1),
      annotation_name_gp = gpar(fontsize = 10),
      show_legend = TRUE
    )

    # Generar Heatmap
    HeatmapEPIC <- Heatmap(combined_data,
                           column_title = "Heatmap EPIC",
                           cluster_rows = TRUE,
                           cluster_columns = TRUE,
                           column_names_gp = gpar(fontsize = 6),
                           row_names_gp = gpar(fontsize = 8),
                           rect_gp = gpar(col = "grey60", lwd = 1),
                           name = 'Z-score',
                           top_annotation = col_ann)

    print(HeatmapEPIC)
  }
  ###### qti
  combined_data <- h_imm_qti
  # Verificar si hay NaN en combined_data
  if (any(is.na(combined_data))) {
    cat("Cannot generate the qti heatmap because there are NaN values in the data.\n")
  } else {
    # Continuar con la generación del heatmap
    combined_data <- h_imm_qti
    col_ann_data <- AnnotData[colnames(combined_data), , drop = FALSE]

    # Asegúrate de que col_ann_data tenga las mismas dimensiones que combined_data
    if (ncol(combined_data) == nrow(col_ann_data)) {
      col_ann_data <- as.matrix(col_ann_data)  # Convertir a matriz
    } else {
      stop("Dimensions of col_ann_data and combined_data do not match")
    }

    # Anotaciones
    col_ann <- HeatmapAnnotation(
      df = col_ann_data[, design_formula, drop = FALSE],  # Selecciona solo la columna design_formula
      gp = gpar(col = "grey60", lwd = 1),
      annotation_name_gp = gpar(fontsize = 10),
      show_legend = TRUE
    )

    # Generar Heatmap
    Heatmap_qti <- Heatmap(combined_data,
                           column_title = "Heatmap qti",
                           cluster_rows = TRUE,
                           cluster_columns = TRUE,
                           column_names_gp = gpar(fontsize = 6),
                           row_names_gp = gpar(fontsize = 8),
                           rect_gp = gpar(col = "grey60", lwd = 1),
                           name = 'Z-score',
                           top_annotation = col_ann)
    print(Heatmap_qti)
  }
  ####### h_imm_xcell
  combined_data <- h_imm_xcell
  # Verificar si hay NaN en combined_data
  if (any(is.na(combined_data))) {
    cat("Cannot generate the xcell heatmap because there are NaN values in the data.\n")
  } else {
    # Continuar con la generación del heatmap
    combined_data <- h_imm_xcell
    col_ann_data <- AnnotData[colnames(combined_data), , drop = FALSE]

    # Asegúrate de que col_ann_data tenga las mismas dimensiones que combined_data
    if (ncol(combined_data) == nrow(col_ann_data)) {
      col_ann_data <- as.matrix(col_ann_data)  # Convertir a matriz
    } else {
      stop("Dimensions of col_ann_data and combined_data do not match")
    }

    # Anotaciones
    col_ann <- HeatmapAnnotation(
      df = col_ann_data[, design_formula, drop = FALSE],
      gp = gpar(col = "grey60", lwd = 1),
      annotation_name_gp = gpar(fontsize = 10),
      show_legend = TRUE
    )

    # Generar Heatmap
    Heatmap_xcell <- Heatmap(combined_data,
                             column_title = "Heatmap xcell",
                             cluster_rows = TRUE,
                             cluster_columns = TRUE,
                             column_names_gp = gpar(fontsize = 6),
                             row_names_gp = gpar(fontsize = 8),
                             rect_gp = gpar(col = "grey60", lwd = 1),
                             name = 'Z-score',
                             top_annotation = col_ann)
    print(Heatmap_xcell)
  }

  cat("\033[32mHeatmap of qti, EPIC and xcell will be stored on plots_convolution_heatmap.pdf\033[0m\n")
  pdf("plots_convolution_heatmap.pdf")
  if (exists("Heatmap_qti")) {
    print(Heatmap_qti)
  } else {
    cat("Heatmap_qti object does not exist.\n")
  }
  if (exists("HeatmapEPIC")) {
    print(HeatmapEPIC)
  } else {
    cat("HeatmapEPIC object does not exist.\n")
  }

  if (exists("Heatmapxcell")) {
    print(Heatmapxcell)
  } else {
    cat("Heatmapxcell object does not exist.\n")
  }
  dev.off()

  # Correlación genes \~ tipos celulares
  ## EPIC
  c_imm_epic <- imm_epic
  c_imm_epic <- as.data.frame(t(c_imm_epic))
  c_imm_epic <- head(c_imm_epic, -1)
  rn_cimmepic <- rownames(c_imm_epic)
  rn_cimmepic[2] <- "CAFs"
  cl_cimmepic <- colnames(c_imm_epic)
  c_imm_epic <- apply(c_imm_epic, 2, as.numeric)
  rownames(c_imm_epic) <- rn_cimmepic
  colnames(c_imm_epic) <- cl_cimmepic
  c_imm_epic <- as.data.frame(c_imm_epic)

  ## Quantiseq
  c_imm_qti <- imm_qti
  c_imm_qti <- as.data.frame(t(c_imm_qti))
  c_imm_qti <- head(c_imm_qti, -1)
  rn_cimmqti <- rownames(c_imm_qti)
  rn_cimmqti[7] <- "T cell CD4+"
  rn_cimmqti[9] <- "T cell regulatory"
  cl_cimmqti <- colnames(c_imm_qti)
  c_imm_qti <- apply(c_imm_qti, 2, as.numeric)
  rownames(c_imm_qti) <- rn_cimmqti
  colnames(c_imm_qti) <- cl_cimmqti
  c_imm_qti <- as.data.frame(c_imm_qti)

  ## xCell
  c_imm_xcell <- as.data.frame(imm_xcell)

  # Renombrar nombres de poblaciones celulares antes de transponer
  c_imm_xcell <- c_imm_xcell %>%
    dplyr::rename('Myeloid DC activated' = 'Myeloid dendritic cell activated') %>%
    dplyr::rename('T cell CD4+' = 'T cell CD4+ (non-regulatory)') %>%
    dplyr::rename('Myeloid DC' = 'Myeloid dendritic cell') %>%
    dplyr::rename('CAFs' = 'Cancer associated fibroblast') %>%
    dplyr::rename('Plasmacytoid DC' = 'Plasmacytoid dendritic cell') %>%
    dplyr::rename('T cell regulatory' = 'T cell regulatory (Tregs)')

  # Transponer
  c_imm_xcell <- as.data.frame(t(c_imm_xcell))
  c_imm_xcell <- head(c_imm_xcell, -4)

  # Poblaciones celulares no interesantes
  xcell_row_delete <- c('Common lymphoid progenitor', 'Common myeloid progenitor',
                        'Granulocyte-monocyte progenitor', 'Hematopoietic stem cell')

  c_imm_xcell <- c_imm_xcell %>%
    filter(!row.names(.) %in% xcell_row_delete)  # Eliminar filas

  rn_cimmxcell <- rownames(c_imm_xcell)
  cl_cimmxcell <- colnames(c_imm_xcell)
  c_imm_xcell <- apply(c_imm_xcell, 2, as.numeric)
  rownames(c_imm_xcell) <- rn_cimmxcell
  colnames(c_imm_xcell) <- cl_cimmxcell
  c_imm_xcell <- as.data.frame(c_imm_xcell)

  decounts_correlation <- function(count_data1, deconv_data, gene_list = NULL, cell_type_list = NULL) {
    # Verificar si las columnas son numéricas
    cat("\033[32mAre numeric?\033[0m\n")
    numeric_columns_check <- sapply(deconv_data, is.numeric)
    print(numeric_columns_check)
    print(table(numeric_columns_check))

    # Verificar si hay valores NA
    cat("\033[32mAre any NA?\033[0m\n")
    na_check <- sapply(deconv_data, function(x) any(is.na(x)))
    print(na_check)
    print(table(na_check))

    if (is.null(gene_list)) {
      gene_list <- rownames(count_data1)
    }

    if (is.null(cell_type_list)) {
      cell_type_list <- rownames(deconv_data)
    }

    cat("\033[32mIniciating the correlation...\033[0m\n")
    correlation_df <- data.frame(matrix(NA, nrow = length(gene_list), ncol = length(cell_type_list), dimnames = list(gene_list, cell_type_list)), check.names = FALSE)
    pvalue_df <- data.frame(matrix(NA, nrow = length(gene_list), ncol = length(cell_type_list), dimnames = list(gene_list, cell_type_list)), check.names = FALSE)

    total_genes <- length(gene_list)
    total_cell_types <- length(cell_type_list)
    for (i in seq_along(gene_list)) {
      gene <- gene_list[i]
      start_time_gene <- Sys.time()
      for (j in seq_along(cell_type_list)) {
        cell_type <- cell_type_list[j]
        # Obtener los valores de expresión del gen y la abundancia del tipo celular
        gene_expression <- as.numeric(count_data1[gene, ])
        cell_type_abundance <- as.numeric(deconv_data[cell_type, ])
        # Realizar la prueba de correlación
        correlation_test <- cor.test(gene_expression, cell_type_abundance, method = "spearman")
        # Almacenar los resultados en los dataframes
        correlation_df[gene, cell_type] <- correlation_test$estimate
        pvalue_df[gene, cell_type] <- correlation_test$p.value
      }
      end_time_gene <- Sys.time()
      cat(sprintf("Done with gene %s (%d/%d). Time taken: %f seconds\n", gene, i, total_genes, as.numeric(difftime(end_time_gene, start_time_gene, units = "secs"))))
    }
    # Retornar los resultados como una lista
    return(list(correlation = correlation_df, pvalues = pvalue_df))
  }


  # Ejecución de la función
  if (correlation){
    cat("\033[32mResults of correlation between gene expression data and cell type abundance data. This data will be imported in .csv\033[0m\n")
    cor_epic <- decounts_correlation(tpm_counts, c_imm_epic)
    write.csv(cor_epic, "correlation_epic.csv")
    cor_qti <- decounts_correlation(tpm_counts, c_imm_qti)
    write.csv(cor_qti, "correlation_qti.csv")
    cor_xcell <-decounts_correlation(tpm_counts, c_imm_xcell)
    write.csv(cor_xcell, "correlation_xcell.csv")}
  data_table_list<- list(EPIC= imm_epic, QTI = imm_qti, XCELL= imm_xcell)
  return(data_table_list)
}
