#' Adds ENTREZID column to data frame
#'
#' Will map a column with ENTREZID for each gene based on the rownames keytype
#' (SYMBOL, ENSEMBL) of data frame. ketype is defaulted to SYMBOL.
#'
#'
#' @param df data frame with SYMBOL or ENSEMBL rownames.
#' @param key keytype, default = "SYMBOL"
#'
#' @import clusterProfiler
#'
#' @importFrom stats na.omit
#' @importFrom AnnotationDbi mapIds
#'
#'
#' @return dataframe with appended ENTREZID column
#' @export
#'

get_Hs_entrez <- function(df, key = "SYMBOL") {

  # Creates a new column with the ENTREZID's in deg_genes
  df$ENTREZID <- mapIds(org.Hs.eg.db,
                                      keys       = rownames(df),
                                      column     = "ENTREZID",
                                      keytype    = key,
                                      multiVals  = "first")

  # There could be N/A values.. not decided if we should omit or not. Must
  # discuss with Luan how to approach. For now i will omit the n/a values.
  df <- na.omit(df)

  return(df)
}
