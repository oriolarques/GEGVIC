#' @title module_ic
#'
#' @description Estimates the composition of immune cell composition from RNA-seq
#' raw data. Additionally it calculates an immunophenogram and immunophenoscores
#' for each sample and groups of study.
#'
#' @param counts Data frame that contains gene expression data as raw counts.
#' @param genes_id Name of the column that contains gene identifiers. Should be
#' one of the following:'entrez_gene_id', 'ensemblgene_id' or 'hgnc_symbol'.
#' @param biomart Data frame containing a biomaRt query with the following
#' attributes: ensembl_gene_id, hgnc_symbol, entrezgene_id, transcript_length,
#' refseq_mrna. In the case of mus musculus data, external_gene_name must be
#' obtained and then change the column name for hgnc_symbol. Uploaded biomaRt
#' queries in GEGVIC: 'ensembl_biomartGRCh37', ensembl_biomartGRCh38_p13' and
#' 'ensembl_biomartGRCm38_p6', 'ensembl_biomartGRCm39'.
#' @param indications Character vector of cancer type codes for each sample in
#' the tpm matrix.This is used by TIMER method. Indications supported
#' can be checked using immunedeconv::timer_available_cancers. Default value is
#' NULL.
#' @param cibersort Path to the CIBERSORT.R and LM22.txt files. Default value is
#' NULL.
#' @param tumor Logical value to define if samples are tumors. If so EPIC and
#' quanTIseq use a signature matrix/procedure optimized for tumor samples.
#' Default value is TRUE.
#' @param rmgenes A character vector of gene symbols. Exclude these genes from
#' the analysis. Use this to exclude e.g. noisy genes.
#' @param scale_mrna Logical. If FALSE, disable correction for mRNA content of
#' different cell types. This is supported by methods that compute an absolute
#' score (EPIC and quanTIseq). Default value is TRUE.
#' @param expected_cell_types Limit the analysis to the cell types given in this
#' list. If the cell types present in the sample are known a priori, setting
#' this can improve results for xCell
#' (see https://github.com/grst/immunedeconv/issues/1).
#' @param metadata Data frame that contains supporting variables to the data.
#' @param response Unquoted name of the variable indicating the groups to analyse.
#' @param compare A character string indicating which method to be used for
#' comparing means. Options are 't.test' and 'wilcox.test' for two groups or
#' 'anova' and 'kruskal.test' for more groups. Default value is NULL.
#' @param p_label Character string specifying label type. Allowed values include
#' 'p.signif' (shows the significance levels), 'p.format' (shows the formatted
#' p-value).
#' @param colors Character vector indicating the colors of the different groups
#' to compare. Default values are two: black and orange.
#'
#' @return Returns ggplot objects showing predicted immune cell populations to
#' be compared between or within samples. Also it returns a list of tables with
#' the data necessary to produce the plots.
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
#' tables_module_ic <- module_ic(counts = input_ge_module,
#'                               genes_id = 'entrezgene_id',
#'                               biomart = ensembl_biomart_GRCh38_p13,
#'                               indications = rep('skcm', ncol(input_ge_module[-1])),
#'                               cibersort = NULL,
#'                               metadata = metadata_ge_module,
#'                               response = Response,
#'                               compare = 'wilcox.test',
#'                               p_label = 'p.format',
#'                               colors = c('black', 'orange'))
#'
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
                      response,
                      compare = NULL,
                      p_label = 'p.format',
                      colors = c('black', 'orange')) {

    # Obtain TPM from raw counts
    print('Obtain TPM from raw counts')

    tpm <- ic_raw_to_tpm(counts = counts,
                         genes_id = genes_id,
                         biomart = biomart)


    # Estimate immune cell fractions
    print('Estimate immune cell fractions')

    ic.pred <- ic_deconv(gene_expression = tpm,
                         indications = indications,
                         cibersort = cibersort,
                         tumor = tumor,
                         rmgenes = rmgenes,
                         scale_mrna = scale_mrna,
                         expected_cell_types = expected_cell_types)


    # Quote the argument response to use it in the next functions
    response <- rlang::enquo(response)

    # Plot a graph comparing immune cell populations between samples group
    print('Estimate immune cell composition')

    ic_plot_comp_samples(df = ic.pred,
                         metadata = metadata,
                         response = !!response,
                         compare = compare,
                         p_label = p_label,
                         colors = colors)

    # Plot a graph comparing immune cell populations within samples
    ic_plot_comp_celltypes(df = ic.pred,
                          metadata = metadata,
                          response = !!response)

    # Calculate and plot immunophenogram (IPG) and immunophenoscores
    ## (IPS) for each sample and each group of study.
    print('Calculate IPG and IPS')

    ips <-ic_score(tpm = tpm,
                   metadata = metadata,
                   response = !!response,
                   compare = compare,
                   p_label = p_label,
                   colors = colors)

    # Return tables generated as a list
    tables_module_ic <- list(tpm = tpm,
                             ic.pred = ic.pred,
                             ips = ips)

    return(tables_module_ic)

}

