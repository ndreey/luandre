#fetches data from a local files:
fetch_data_from_file <- function(counts, sample_table) {
  #declares count data with a global variable
  count_data <<- read.table(
    file = counts,
    header = TRUE,
    row.names = 1,
    as.is = TRUE,
    sep = "\t")

  #declares sample table with a global variable:
  sample_table <<- read.table(
    file = sample_table,
    header = TRUE,
    row.names = 1,
    as.is = TRUE,
    sep = "\t")
}

