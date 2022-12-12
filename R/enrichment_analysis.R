enrich_analysis <- function(counts_file, sample_table_file, db, ontology = NULL){
  fetch_data_from_file(counts_file, sample_table_file)
  digital_gene_expression(count_data, sample_table)
  deg <- fit_test_model(dge_list, design_matrix)
  filtered_degs <- get_deg(deg)
  hs_entrez <- get_Hs_entrez(filtered_degs)
  res <- enrich_Hs_genes(hs_entrez$ENTREZID, db, ontology)
  enrich_plots(res, 10, font_size = 10, file_id = "EA",width = 8, height = 6, dpi = 300)
}
