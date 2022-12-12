
#' Title
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


fetch_data_from_file("E-MTAB-2523.counts.txt", "E-MTAB-2523_sample table.txt")

digital_gene_expression <- function(count_data, sample_table) {

  #Filters count data:
  meanlog2cpm <- rowMeans(log2(cpm(count_data)+1))
  count_data <- count_data[meanlog2cpm > 1,]

  #factorize control/disease. Requires column with the name "disease"
  control_factor <- factor(sample_table$disease)
  designDF <- data.frame(controls = control_factor)
  design_matrix <<- model.matrix( ~ 0 + controls , data = designDF)
  colnames <- colnames(design_matrix)

  dge <- DGEList(count_data)
  dge_list <<- calcNormFactors(dge)

  case <<- colnames[1]
  control <<- colnames[2]
}

