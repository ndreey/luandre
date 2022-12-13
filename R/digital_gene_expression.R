digital_gene_expression <- function(count_data, sample_table) {

  #-----------------------------------------------------------------------------
  # CAN WE HAVE A MORE VIGORUS FILTERING?
  #
  # Normalizes all the counts with counts per million (CPM)
  # --> cpm(count_data)                       # makes it cpm
  #
  # Sums how many of the cpm values for each gene and sample that were > 100.
  # --> rowSums(cpm(count_data) > 100)        # sums how many cpm > 100 per gene
  #
  # Keeps the genes that at least have two samples with cpm > 100 by indexing.
  # keep <- rowSums(cpm(count_data>100) >= 2  # index of TRUE or FALSE
  # counts <- count_data[keep,]
  #
  # This will give us 3834 possible DEGs instead of 13888.
  #-----------------------------------------------------------------------------
  #Filters count data:
  meanlog2cpm <- rowMeans(log2(cpm(count_data)+1))
  count_data <- count_data[meanlog2cpm > 1,]
  # print(dim(count_data))     # 13888


  #-----------------------------------------------------------------------------
  # CHECK LINE 31
  # Missed setting levels = c('normal', 'carcinoma') in factorizing.Important
  # that we have control as 0 (base line) so we can see how much
  # disease diverges.
  #-----------------------------------------------------------------------------
  #factorize control/disease. Requires column with the name "disease"
  control_factor <- factor(sample_table$disease)
  designDF <- data.frame(controls = control_factor)


  #-----------------------------------------------------------------------------
  # CHECK LINE 42
  # I believe you meant to do colnames(levels(control_factor)) to set colnames
  # to carcinoma and normal.
  # ----------------------------------------------------------------------------
  # Creates design_matrix
  design_matrix <<- model.matrix( ~ 0 + controls , data = designDF)
  colnames <- colnames(design_matrix)


  #-----------------------------------------------------------------------------
  # CHECK LINE 50
  # Because group = control_factor was not added, dge$samples only consists of
  # 1's. When we add the grouping we also get groups into this DGEList object.
  #-----------------------------------------------------------------------------
  dge <- DGEList(count_data)
  dge_list <<- calcNormFactors(dge)


  case <<- colnames[1]
  control <<- colnames[2]
}
