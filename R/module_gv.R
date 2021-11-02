#' @titlte module_gv
#'
#' @description
#'
#' @param muts
#' @param metadata
#' @param response
#' @param top_genes
#' @param specific_genes
#' @param colors
#' @param compare
#' @param p_label
#' @param gbuild
#' @param mut_sigs
#' @param tri.counts.method
#'
#' @return
#'
#' @export
#'
#' @importFrom maftools read.maf
#' @importFrom maftools plotmafSummary
#' @importFrom maftools oncoplot
#' @import rlang
#' @import dplyr
#' @import ggplot2
#' @import ggpubr
#' @import tibble
#' @import tidyr
#' @import deconstructSigs
#' @import pheatmap
#' @import ggplotify
#'
#' @examples
#' module_gv(muts = input_gv_module,
#'           metadata = metadata_ge_module,
#'           response = Response,
#'           top_genes = 10,
#'           specific_genes = NULL,
#'           colors = c('black' ,'orange'),
#'           compare = 'wilcox.test',
#'           p_label = 'p.format',
#'           gbuild = 'BSgenome.Hsapiens.UCSC.hg19',
#'           mut_sigs = COSMIC_v2_SBS_GRCh37,
#'           tri.counts.method = 'default')
#'
module_gv <- function(muts,
                      metadata,
                      response,
                      top_genes = 10,
                      specific_genes = NULL,
                      colors = c('black' ,'orange'),
                      compare = NULL,
                      p_label = 'p.format',
                      gbuild = 'BSgenome.Hsapiens.UCSC.hg19',
                      mut_sigs = COSMIC_v2_SBS_GRCh37,
                      tri.counts.method = 'default') {

    response <- enquo(response)

    # Show genetic variations summary
    gv_mut_summary(muts = muts,
                   metadata = metadata,
                   response = !!response,
                   top_genes = top_genes,
                   specific_genes = specific_genes,
                   colors = colors)

    # Calculate mutational load
    gv_mut_load(muts = muts,
                metadata = metadata,
                response = !!response,
                compare = compare,
                p_label = p_label,
                colors = colors)

    # Extract mutational signature profiles
    gv_mut_signatures(muts = muts,
                      metadata = metadata,
                      response = !!response,
                      gbuild = gbuild,
                      mut_sigs = mut_sigs,
                      tri.counts.method = tri.counts.method,
                      colors = colors)



}
