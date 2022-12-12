
#' Title Fetch data from file
#'
#' @param counts
#' @param sample_table
#'
#' @return
#' @export
#'
#' @import edgeR
#'
#' @examples
#fetches data from a local files:
fetch_data_from_file <- function(counts_file, sample_table_file) {

  #declares count data with a global variable:
  count_data <<- read.table(
    file = counts_file,
    header = TRUE,
    row.names = 1,
    as.is = TRUE,
    sep = "\t")

  #declares sample table with a global variable:
  sample_table <<- read.table(
    file = sample_table_file,
    header = TRUE,
    row.names = 1,
    as.is = TRUE,
    sep = "\t")
}






