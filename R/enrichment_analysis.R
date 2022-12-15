#'Enrich Analysis.
#'
#' Calls upon all functions of the package to instantly perform EA and save the result plots to the users working directory.
#'
#' @param counts_file file containing count data.
#' @param sample_table_file file containing sample table.
#' @param db "KEGG" or "GO"
#' @param ontology "BP", "CC" or "MF"
#' @param threshold minimum mean-log2-CPM to keep the data during low expression filtering.
#'
#' @export


enrich_analysis <- function(counts_file, sample_table_file, db, ontology = NULL, threshold = 1){
  data <- fetch_data_from_file(counts_file, sample_table_file)
  data <- filter(data, threshold)
  dge <- digital_gene_expression(data)
  deg <- fit_test_model(dge[[2]], dge[[1]])
  filtered_degs <- get_deg(deg)
  hs_entrez <- get_Hs_entrez(filtered_degs)
  rich_res <- enrich_Hs_genes(hs_entrez$ENTREZID, db, ontology)
  enrich_plots(rich_res, 10, font_size = 10, file_id = "EA",width = 8, height = 6, dpi = 300)
  save_deg_excel(filtered_degs, "DEGs")
}



