#' HTG_calculate_summary_stats
#'
#' @description This function calculates summary statistics for each column of the input data frame, including minimum, maximum, mean, median, mode, standard deviation, variance, range, first quartile (Q1), third quartile (Q3), interquartile range (IQR), skewness, kurtosis, count of missing values, and coefficient of variation (CV).
#'
#' @param counts_data A data frame containing the input counts.
#' @param pattern (Optional) A regular expression pattern to identify control probes in the count data. For HTG, this could be "^NC-|^POS-|^GDNA-|^ERCC-". If NULL, the pattern will not be applied.
#'
#' @return A data frame containing summary statistics for each column of the input data. and also will be .csv
#' @export
#'
#' @importFrom utils write.csv
#'
#' @examples
#' summary <- HTG_calculate_summary_stats(counts_data_tutorial, pattern = "^NC-|^POS-|^GDNA-|^ERCC-")
#' @name HTG_calculate_summary_stats
HTG_calculate_summary_stats <- function(counts_data, pattern = NULL) {
  if (!is.null(pattern)) {
    counts_filtered <- subset(counts_data, !grepl(pattern, rownames(counts_data)))
  }
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
  write.csv(summary_stats, file = "summary_stats.csv")
  return(summary_stats)
}
