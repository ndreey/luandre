BI503G R-Package Build Seminar
================
André Bourbonnais, Luan Gardenalli
2022-12-15

<!-- README.md is generated from README.Rmd. Please edit that file -->

# LuAndre

<!-- badges: start -->
<!-- badges: end -->

Package created for Differential Expression Analysis of bulk RNA-seq
data.

# Installation

LuAndre is not available on CRAN: Install the development version of
LuAndre from [GitHub](https://github.com/ndreey/luandre) use:

``` r
# install.packages("devtools")
devtools::install_github("ndreey/luandre")
```

Manual installment of Bioconductor packages might be needed.

    if(BiocManager::texttocome)

## Tutorials

- Check out the [LuAndre Manual](https://raw.githack.com/ndreey/seminar/main/docs/index.html) for more information

## Enrichment Analysis

### Option 1: Use the wrapper

This wrapper orchestrates the necessary functions calls in order

``` r
# Run Harmony to integrate the reference cells
enrich_analysis("E-MTAB-2523.counts.txt",   # counts table
                "E-MTAB-2523.sample.txt",   # sample table
                db = "GO",                  # "KEGG" or "GO"
                ontology = "BP",            # Biological Process
                threshold = 2)              # Filtering threshold
```

Note that `ontology =` is only required if `"GO"`is set.

### Option 2: Step-wise analysis

The plots and excel file will be stored in your working directory.

``` r
# Import data
fetched_data <- fetch_data_from_file("E-MTAB-2523.counts.txt",
                                   "E-MTAB-2523.sample.txt")

# Lets check the dims
dim(fetched_data[[1]])       # 19351 genes, 18 samples
dim(fetched_data[[2]])       # 18 samples, 3 sample info columns.

# Lets filter the data to remove low expression counts
filtered_data <- filter(fetched_data, threshold = 3)

# Lets check how many gene are kept
dim(filtered_data[[1]])       # 10730 genes

# Lets create a Digital gene expression list and a design matrix using data3
digital_list <- digital_gene_expression(filtered_data)

# Lets fit it to a general linear model
results <- fit_test_model(digital_list)

# Lets filter out degs that dont meet our requirments.
deg_df <- get_deg(results, lfc = 2, adj_p = 0.01)
dim(deg_df)      # 180 DEGs!

# Lets store them to an excel file.
save_deg_excel(deg_df3, "E-MTAB-2523")

# Prepare for plotting
# Adds ENTREZID column to data frame
deg_df <- get_Hs_entrez(deg_df)


# Enrich genes with KEGG or GO (Homo sapiens)
enrich_results <- enrich_Hs_genes(deg_df$ENTREZID, db = "GO", ontology = "BP")

# Plots for over-representaiton data.
enrich_plots(enrich_results, n = 10, font_size = 10, file_id = "E-MTAB-2523", 
             width = 8, height = 6, dpi = 300)
```

# Installation notes

## System requirements:

LuAndre has been successfully installed on windows 11 using the devtools
package to install from GitHub.

Dependencies:

- R\>=3.6.x
- edgeR
- limma
- clusterProfiler
- org.Hs.eg.db
- ggplot2
- enrichplot
- GOSemSim
- ggupset
- openxlsx
- AnnotationDbi

## Troubleshooting:

- You may need to install the latest version of devtools (because of the
  recent GitHub change from “master” to “main” terminology, which can
  cause `install_github` to fail).
