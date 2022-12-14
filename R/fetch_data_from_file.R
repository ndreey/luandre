#'Fetch data from file.
#'
#' Imports transcriptomics count data from local text files.
#' Expects two files, a count data file and a sample table file.
#'
#' @param counts_file file containing count data.
#' @param sample_table_file file containing sample table.
#'
#' @importFrom utils read.table
#'
#' @return list with a count data-frame and a sample table data-frame.
#' @export

#fetches data from a local files:
fetch_data_from_file <- function(counts_file, sample_table_file) {

  #declares count data
  count_data <- read.table(
    file = counts_file,
    header = TRUE,
    row.names = 1,
    as.is = TRUE,
    sep = "\t")

  #declares sample table
  sample_table <- read.table(
    file = sample_table_file,
    header = TRUE,
    row.names = 1,
    as.is = TRUE,
    sep = "\t")

  return(list(count_data, sample_table))
}





