#' Plots dot-, bar-, cnet-, upset-, go- and treeplot for
#' over-representaiton data.
#'
#' Each plot will be saved to a .png following this naming convention:
#' (file id)_(database type)_(plot type).png.
#' Database type and plot type is set automatically, user must set unique id.
#' User also has the ability to change the font size and number of annotations.
#' Furthermore, user can decide width, depth and dpi of .png file by passing the
#' dots. The plots are saved in the current working directory.
#'
#'
#' @param rich_res enrichResult class
#' @param n number of annotations
#' @param font_size font size
#' @param file_id unique id for plots
#' @param ... width =, height =, dpi =,
#'
#' @import enrichplot
#' @import GOSemSim
#' @import ggplot2
#'
#' @return null
#' @export
#'
enrich_plots <- function(rich_res, n, font_size = NULL, file_id, ... ) {

  # Idea------------------------------------------------------------------------
  # Would be real cool if we could split plot so you have control vs case.
  # But we can do that in later versions.Lets just get basic plotting working.


  # Store the ontology of enrichRes object
  type <- toupper(rich_res@ontology)

  # Define a string with placeholders for multiple values
  plot_title <- "%s - %s"

  # Replace the placeholders with real values
  if (type == "KEGG"){

    plot_title <- sprintf(plot_title, "KEGG", "Enriched Pathways")

  } else if (type == "BP") {

    plot_title <- sprintf(plot_title, "GO", "Biological Processes")

  } else if (type == "CC") {

    plot_title <- sprintf(plot_title, "GO", "Cellular Components")

  } else if (type == "MF") {

    plot_title <- sprintf(plot_title, "GO", "Molecular Functions")

  }

  # Draw plots------------------------------------------------------------------
  #

  # Dotplot
  dot <- enrichplot::dotplot(rich_res, showCategory = n, title = plot_title,
                             font.size = font_size)

  # Barplot
  bar <- enrichplot::barplot(rich_res, showCategory = n, title = plot_title,
                             font.size = font_size)

  # Network plot
  cnet <- enrichplot::cnetplot(rich_res, colorEdge = TRUE) +
    ggtitle(ggtitle("Gene-Concept Network", subtitle = plot_title))

  # Not 100% sure what this is...
  upset <- enrichplot::upsetplot(rich_res) +
    ggtitle(ggtitle("Upset plot", subtitle = plot_title))

  # Like CNET but we get arrows!
  goplot <- enrichplot::goplot(rich_res) +
    ggtitle(ggtitle("GOplot", subtitle = plot_title))


  # Nested list holding plots generated and the name.
  plots <- list(list("dot",dot), list("bar", bar), list("cnet", cnet),
                list("upset", upset))

  # Treeplot--------------------------------------------------------------------
  #

  # Creates GOSemSimDATA
  if (type != "KEGG"){

    semsim <- GOSemSim::godata("org.Hs.eg.db", ont = type)
    sim_mat <- enrichplot::pairwise_termsim(rich_res, semData = semsim)

    # Plots treeplot
    tree <- enrichplot::treeplot(sim_mat, showCategory = n) +
      ggtitle(ggtitle("Treeplot", subtitle = plot_title))

    # Adds the tree plot to list
    plots <- append(plots, list("tree", tree))

  } else {
    # Do not plot tree for KEGG.
  }

  # Plot to png-----------------------------------------------------------------
  # For each plot created, we save them as .png with unique identifiers.
  #

  # For each plot we save them as (unique id)_(database type)_(plot type).png
  for (plot in plots) {
    file_name <- "%s_%splot_%s.png"
    file_name <- sprintf(file_name, file_id, plot[1], type)
    print(file_name)

    # User can "pass in the dots" of width, height and dpi.
    ggplot2::ggsave(file_name, plot[[2]], bg = "white", ...)
  }

}