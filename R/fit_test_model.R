
#' Linear regression using generalized linear model
#'
#' Estimates common, trended and tagwise dispersions and fits a general
#' linear model to test for differentially expressed genes using functions from
#' the edgeR and limma package
#'
#' @param dge_list DGEList object
#' @param design_matrix numeric design matrix
#'
#' @import edgeR
#' @import limma
#'
#' @return Digital Gene Expression Likelihood Ratio Test data and results
#' - class (DGELRT)
#' @export
#'

fit_test_model <- function(dge_list, design_matrix) {

  # Estimate Common, Trended and Tagwise dispersions and adds to dge_list
  dge_list <- edgeR::estimateDisp(dge_list, design_matrix)

  # Fit the model
  fit <- edgeR::glmQLFit(dge_list, design_matrix)

  # Contrast matrix
  contrast_matrix <- limma::makeContrasts(CaseVsControl = case - control,
                                levels = design_matrix)

  # Statistical testing
  deg_results <- edgeR::glmQLFTest(fit, contrast = contrast_matrix)

  # Change name "-1*samp_groupnormal 1*samp_groupcase" to "Control vs Case"
  deg_results$comparison <- "Control vs Case"

  # summary without lfc or pvalue limit
  # Prints out Up and Down
  summary(limma::decideTests(deg_results))

  return(deg_results)
}


