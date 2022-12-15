#'Enrich Analysis.
#'
#' Calls upon all functions of the package to instantly perform EA and save the
#' result plots to the users working directory.
#'
#' @param counts_file file containing count data.
#' @param sample_table_file file containing sample table.
#' @param db "KEGG" or "GO"
#' @param ontology "BP", "CC" or "MF"
#' @param threshold minimum mean-log2-CPM to keep the data during low
#' expression filtering.
#'
#' @export


enrich_analysis <- function(counts_file, sample_table_file, db, ontology = NULL,
                            threshold = 1){
  # Imports transcriptomics count data from local text files.
  data <- fetch_data_from_file(counts_file, sample_table_file)

  # Filter low expressing counts by mean logarithm counts-per-million.
  data <- filter(data, threshold)

  # Creates a Digital gene expression list and a design matrix.
  dge <- digital_gene_expression(data)

  # Linear regression using generalized linear model
  deg <- fit_test_model(data)

  # Returns a data frame with DEGS meeting the inputed log fold change and
  # adjusted p-values (FDR) conditions.
  filtered_degs <- get_deg(deg)

  # Adds ENTREZID column to data frame
  hs_entrez <- get_Hs_entrez(filtered_degs)

  # Enrich genes with KEGG or GO (Homo sapiens)
  rich_res <- enrich_Hs_genes(hs_entrez$ENTREZID, db, ontology)

  # Plots for over-representaiton data.
  enrich_plots(rich_res, 10, font_size = 10, file_id = "EA",width = 8,
               height = 6, dpi = 300)

  # Save the filtered_degs data.frame  to excel file.
  save_deg_excel(filtered_degs, "DEGs")
}



