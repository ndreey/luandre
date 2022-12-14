#' Enrich genes with KEGG or GO pathways (Homo sapiens only)
#'
#' Returns a enrichResult class using a vector of ENTREZID's and specifying
#' which database to enrich from (KEGG or GO) using enrichGO or enrichKEGG
#' function from clusterProfiler. If db = "GO", ontology must be chosen.
#' Biological Processes ("BP"), Cellular Components ("CC") or Molecular
#' Functions ("MF").
#'
#' @param entrez vector with ENTREZID's to enrich
#' @param db "KEGG" or "GO"
#' @param ontology "BP", "CC" or "MF"
#'
#' @import clusterProfiler
#' @import org.Hs.eg.db
#'
#' @return enrichResult class object
#' @export
#'
enrich_Hs_genes <- function(entrez, db, ontology = NULL) {

  # Print statement
  print("Loading annotations...")

  if (db == "GO") {
    # Enrichment GO
    rich_res <- enrichGO(gene          = entrez,
                         OrgDb         = org.Hs.eg.db,
                         ont           = ontology,
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05)

  } else if (db == "KEGG") {
    # Enrichment KEGG
    rich_res <- enrichKEGG(gene          = entrez,
                           organism      = "hsa",
                           keyType       = "kegg",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.05)
  }

  # print stamtent letting user know it is complete
  print("Annotations loaded")

  return(rich_res)
}
