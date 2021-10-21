#' @title module_ic
#'
#' @param counts
#' @param genes_id
#' @param biomart
#' @param indications
#' @param cibersort
#' @param tumor
#' @param rmgenes
#' @param scale_mrna
#' @param expected_cell_types
#' @param metadata
#' @param response
#' @param compare
#' @param p_label
#' @param colors
#'
#' @return
#'
#' @export
#'
#' @import dplyr
#' @import tibble
#' @import immunedeconv
#' @import tidyr
#' @import ggplot2
#' @import ggpubr
#' @import rlang
#'
#' @examples
module_ic <- function(counts,
                      genes_id,
                      biomart,
                      indications = NULL,
                      cibersort = NULL,
                      tumor = TRUE,
                      rmgenes = NULL,
                      scale_mrna = TRUE,
                      expected_cell_types = NULL,
                      metadata,
                      response, # unquoted name of grouping variable
                      compare = NULL, # compare method ggpubr: t.test, wilcox.test, anova, kruskal.test
                      p_label = 'p.format',
                      colors = c('black', 'orange')) {

    # Obtain TPM from raw counts
    tpm <- ic_raw_to_tpm(counts = counts,
                         genes_id = genes_id,
                         biomart = biomart)

    # Estimate immune cell fractions
    ic.pred <- ic_deconv(gene_expression = tpm,
                         indications = indications,
                         cibersort = cibersort,
                         tumor = tumor,
                         rmgenes = rmgenes,
                         scale_mrna = scale_mrna,
                         expected_cell_types = expected_cell_types)


    # Quote the argument response to use it in the next functions
    response <- enquo(response)

    # Plot a graph comparing immune cell populations between samples group
    ic_plot_comp_samples(df = ic.pred,
                         metadata = metadata,
                         response = !!response, # unquoted name of grouping variable
                         compare = compare, # compare method ggpubr: t.test, wilcox.test, anova, kruskal.test
                         p_label = p_label,
                         colors = colors)

    # Plot a graph comparing immune cell populations within samples
    ic_plot_comp_celltypes(df = ic.pred,
                          metadata = metadata,
                          response = !!response)

}

