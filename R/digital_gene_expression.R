#'Digital Gene Expression.
#'
#' Creates a Digital gene expression list and a design matrix.
#'
#' @param data list consisting of the count data and sample table.
#'
#' @importFrom stats model.matrix
#'
#' @return list containing a design matrix and a DGElist object.
#' @export

digital_gene_expression <- function(data) {

  count_data <- data[[1]]
  sample_table <- data[[2]]

  #factorize control/disease. Requires column with the name "disease"
  control_factor <- factor(sample_table$disease)
  designDF <- data.frame(controls = control_factor)
  design_matrix <- model.matrix( ~ 0 + controls , data = designDF)

  dge <- DGEList(count_data)
  dge_list <- calcNormFactors(dge)

  colnames(design_matrix)[1] = "cases"
  colnames(design_matrix)[2] = "control"

  return(list(design_matrix, dge_list))
}
