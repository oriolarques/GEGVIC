#' @title ic_deconv
#'
#' @description Estimates immune cell fractions in RNA-seq expression data using
#' the following methods QUANTISEQ, TIMER, MCP_COUNTER, XCELL, EPIC and
#' CIBERSORT (absolute mode). To use the latter, the user need to register on
#' the CIBERSORT web page (https://cibersort.stanford.edu), obtain a license and
#' download the source code in form of two files CIBERSORT.R and LM22.txt. Then,
#' the user need to specify the path to the storage location of such files in
#' the cibersort parameter.
#'
#' @param gene_expression Output from the ic_raw_to_tpm function. This is a
#' matrix containing expression counts as TPM with HGNC gene symbols as rownames
#' and samples identifiers as colnames.
#' @param indications Vector used by TIMER method that specifies a cancer
#' indication for each sample in the tpm matrix. Indications supported by TIMER
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
#'
#' @return Returns a data frame.
#'
#' @export
#'
#' @import immunedeconv
#' @import dplyr
#' @import tidyr
#'
#' @examples
#' tpm <- ic_raw_to_tpm(counts = input_ge_module,
#'                      genes_id = 'entrezgene_id,
#'                      biomart = ensembl_biomart_GRCh38_p13)
#' ic.pred <- ic_deconv(gene_expression = tpm,
#'                      indications = rep('skcm', ncol(tpm)),
#'                      cibersort = 'cibersort/',
#'                      tumor = TRUE,
#'                      rmgenes = NULL,
#'                      scale_mrna = TRUE,
#'                      expected_cell_types = NULL)

ic_deconv <- function(gene_expression,
                      indications = NULL,
                      cibersort = NULL,
                      tumor = TRUE,
                      rmgenes = NULL,
                      scale_mrna = TRUE,
                      expected_cell_types = NULL) {

    library(immunedeconv)  # load immunedeconv to avoid Error: object 'xCell.data' not found


    # Define immunedeconv methods except CIBERSORT
    idc_methods <- c('quantiseq', 'timer',
                     'mcp_counter', 'xcell', 'epic')

    # Create necessary objects
    temp_df <- NULL # Data frame to store the predictions for each method
    ic_pred <- NULL # Final data frame to store all predictions

    # Iterate over the methods vector
    for (i in seq_along(idc_methods)) {
        # Obtain predictions for each method
        temp_df <- immunedeconv::deconvolute(gene_expression = gene_expression,
                                             method = idc_methods[i],
                                             indications = indications,
                                             tumor = tumor)
        # Add the method name to the cell_type column
        temp_df <- temp_df %>%
            dplyr::mutate(cell_type = paste0(cell_type, '_', idc_methods[i]))

        # Store predictions in ic_pred object
        ## In the first iteration ic_pred = temp_df
        if (i == 1) {
            ic_pred <- temp_df
        ## Afterwards add rows to ic_pred object
        } else {
            ic_pred <- bind_rows(ic_pred, temp_df)
        }

    }

    # If the user has access to CIBERSOFT
    if (is.null(cibersort) == FALSE){
        # Set the path to CIBERSORT and its matrix
        set_cibersort_binary(paste0(cibersort, 'CIBERSORT.R'))
        set_cibersort_mat(paste0(cibersort, 'LM22.txt'))

        # Obtain predictions for CIBERSORT
        temp_df <- immunedeconv::deconvolute(gene_expression = gene_expression,
                                             method = 'cibersort')
        # Add the method name to the cell_type column
        temp_df <- temp_df %>%
            dplyr::mutate(cell_type = paste0(cell_type, '_CIBERSORT'))

        # Add predictions in ic_pred object
        ic_pred <- bind_rows(ic_pred, temp_df)

        # Obtain predictions for CIBERSORT-ABS
        temp_df <- immunedeconv::deconvolute(gene_expression = gene_expression,
                                             method = 'cibersort_abs')
        # Add the method name to the cell_type column
        temp_df <- temp_df %>%
            dplyr::mutate(cell_type = paste0(cell_type, '_CIBERSORT-ABS'))

        # Add predictions in ic_pred object
        ic_pred <- bind_rows(ic_pred, temp_df)

    }

    # Separate method name from cell_type
    ic_pred <- ic_pred %>%
        tidyr::separate(cell_type, into = c('cell_type', 'method'), sep = '_') %>%
        dplyr::mutate(method = toupper(method)) %>%
        # Set MCP to MCP_COUNTER
        dplyr::mutate(method = ifelse(method == 'MCP', 'MCP_COUNTER', method))

    return(ic_pred)

}
