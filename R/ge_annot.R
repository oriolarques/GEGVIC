#' @title ge_annot
#'
#' @description
#'
#' @param results_dds The output of the ge_diff_exp function. List containing
#' data frames of differential gene expression results between the different groups.
#' @param genes_id Name of the column that contains gene identifiers. Should be
#' one of the following: 'ensembl_gene_id', 'entrezgene_id' or 'hgnc_symbol'.
#'
#' @return Returns a list of data frames containing differential gene expression
#' results, one for each level comparison where genes have been annotated using
#' the hgnc_symbol: Hugo Gene Nomenclature Committee nomenclature.
#'
#' @export
#'
#' @import DESeq2
#' @import dplyr
#' @import tibble
#' @import rlang
#'
#' @examples
#' ge_annot(results_dds = results_dds, genes_id = 'entrezgene_id')

ge_annot <- function(results_dds,
                     genes_id){

    # Create a data.frame to store temporarilly data
    temp_df <- NULL

    # Create a list to store the annotated results
    results_dds.annot <- list()

    # Iterate over the differential gene expression results list
    for (i in seq_along(results_dds)) {

        genes_id <- enquo(genes_id)

        # Get the results as a data frame
        # First evaluate if genes are annotated already as hgnc_symbol
        if (rlang::as_name(genes_id) == 'hgnc_symbol'){
            # If so:
            temp_df <- as.data.frame(results_dds[[i]]) %>%
                # Rownames to column with the name indicated in the genes_id parameter
                tibble::rownames_to_column(rlang::as_name(genes_id)) %>%
                # Filter those missing gene symbols
                dplyr::filter(hgnc_symbol != '') %>%
                # Remove duplicated genes
                dplyr::distinct(hgnc_symbol, .keep_all = TRUE)

        } else {
            # If genes are identified as entrezgene_id or ensembl_gene_id:
            temp_df <- as.data.frame(results_dds[[i]]) %>%
                # Rownames to column with the name indicated in the genes_id parameter
                tibble::rownames_to_column(rlang::as_name(genes_id)) %>%
                # Join the data frame with the GRCh38_p13 biomaRt table stored in data
                dplyr::inner_join(x = .,
                           y = ensembl_biomart_GRCh38_p13 %>%
                               mutate(entrezgene_id = as.character(entrezgene_id)),
                           by = rlang::as_name(genes_id)) %>%
                # From all annotation columns keep only hgnc symbol column
                dplyr::select(hgnc_symbol,
                              everything(),
                              -c(entrezgene_id, ensembl_gene_id)) %>%
                # Filter those missing gene symbols
                dplyr::filter(hgnc_symbol != '') %>%
                # Remove duplicated genes
                dplyr::distinct(hgnc_symbol, .keep_all = TRUE)
        }
        # Store the results in the list
        results_dds.annot[[i]] <- temp_df
        names(results_dds.annot)[i] <- names(results_dds)[i]

    }

    return(results_dds.annot)
}
