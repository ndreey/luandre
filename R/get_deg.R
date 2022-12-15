#' Filters out differently expressed genes (DEGs).
#'
#' Returns a data frame with DEGS meeting the inputed log fold change and
#' adjusted p-values (FDR) conditions.
#'
#' @param deg_results DGELRT object
#' @param lfc numeric value indicating min log fold change. Default = 1
#' @param adj_p numeric value indicating max adjusted p value. Default 0.05
#'
#' @import edgeR
#' @import limma
#'
#' @return data.frame with the DEGs
#' @export
#'
#'
get_deg <- function(deg_results, lfc = 1, adj_p = 0.05) {

  # Sort the genes based on p-value
  deg_table <- topTags(deg_results, adjust.method = "BH",
                              n = nrow(deg_results$table))

  # Filter out degs based on user input
  deg_genes <- subset(deg_table$table, FDR < adj_p & abs(logFC) > lfc)

  # Return the filtered table
  return(deg_genes)
}
