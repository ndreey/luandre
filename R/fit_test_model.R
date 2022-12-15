
#' Linear regression using generalized linear model
#'
#' Estimates common, trended and tagwise dispersions and fits a general
#' linear model to test for differentially expressed genes using functions from
#' the edgeR and limma package
#'
#' @param digital_deg_list list containing a design matrix and a DGElist object.
#'
#' @import edgeR
#' @import limma
#'
#' @return Digital Gene Expression Likelihood Ratio Test data and results
#' - class (DGELRT)
#' @export
#'

fit_test_model <- function(digital_deg_list) {

  design_matrix <- digital_deg_list[[1]]
  dge_list <- digital_deg_list[[2]]

  # Estimate Common, Trended and Tagwise dispersions and adds to dge_list
  dge_list <- estimateDisp(dge_list, design_matrix)

  # Fit the model
  fit <- glmQLFit(dge_list, design_matrix)

  # Contrast matrix
  contrast_matrix <- makeContrasts(CaseVsControl = cases - control,
                                levels = design_matrix)

  # Statistical testing
  deg_results <- glmQLFTest(fit, contrast = contrast_matrix)

  # Change name "-1*samp_groupnormal 1*samp_groupcase" to "Control vs Case"
  deg_results$comparison <- "Control vs Case"


  return(deg_results)
}


