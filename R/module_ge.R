#' @title module_ge
#'
#' @description
#'
#' @param counts
#' @param genes_id
#' @param metadata
#' @param colors
#' @param design
#' @param ref_level
#' @param fold_change
#' @param p.adj
#' @param gmt
#' @param gsea_pvalue
#'
#' @return
#'
#' @export
#'
#' @import DESeq2
#' @import dplyr
#' @import ggplot2
#' @import tibble
#' @import apeglm
#' @import ggrepel
#' @import clusterProfiler
#' @import GSEAmining
#'
#' @examples
#'
#'
module_ge <- function(counts,
                      genes_id,
                      metadata,
                      colors = c('black', 'orange'),
                      design,
                      ref_level,
                      biomart,
                      fold_change = 2,
                      p.adj = 0.05,
                      gmt,
                      gsea_pvalue = 0.2) {

    # Get data PCA
    pca <- ge_pca(counts,
           genes_id,
           metadata,
           design,
           colors = c('black', 'orange'))

    print(pca)

    # Run differential gene expression analysis
    results.dds <- ge_diff_exp(counts = counts,
                               genes_id = genes_id,
                               metadata = metadata,
                               design = design,
                               ref_level = ref_level)


    # Annotate gene symbols
    annot.res <- ge_annot(results_dds = results.dds,
                          genes_id = genes_id,
                          biomart = biomart)

    # Create volcano plot
    ge_volcano(annot_res = annot.res,
               fold_change = fold_change,
               p.adj = p.adj)

    # Obtain GSEA results
    ge_gsea(annot_res = annot.res,
            gmt = gmt,
            gsea_pvalue = gsea_pvalue)


}

