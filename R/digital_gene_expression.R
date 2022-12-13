digital_gene_expression <- function(count_data, sample_table) {

  #Filters count data:
  meanlog2cpm <- rowMeans(log2(cpm(count_data)+1))
  count_data <- count_data[meanlog2cpm > 1,]
  # print(dim(count_data))     # 13888

  #factorize control/disease. Requires column with the name "disease"
  control_factor <- factor(sample_table$disease)
  designDF <- data.frame(controls = control_factor)

  # Creates design_matrix
  design_matrix <<- model.matrix( ~ 0 + controls , data = designDF)
  colnames <- colnames(design_matrix)

  dge <- DGEList(count_data)
  dge_list <<- calcNormFactors(dge)


  case <<- colnames[1]
  control <<- colnames[2]
}
