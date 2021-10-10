#' Title
#'
#' @param results_dds The output of the ge_diff_exp function. List containing
#' differential gene expression results between the different groups.
#' @param genes_id Name of the column that contains gene identifiers. Should be
#' one of the following:'entrez_gene_id', 'ensembl_gene_id' or 'hgnc_symbol'.
#'
#' @return
#' @export
#'
#' @import DESeq2
#' @import dplyr
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
    for(i in seq_along(results_dds)){

        genes_id <- enquo(genes_id)

        # Get the results as a data frame
        temp_df <- as.data.frame(results_dds[[i]]) %>%
            # Rownames to column with the name indicated in the genes_id parameter
            rownames_to_column(as_name(genes_id)) %>%
            inner_join(x = .,
                       y = ensembl_biomart_GRCh38_p13 %>%
                           mutate(entrezgene_id = as.character(entrezgene_id)),
                       by = as_name(genes_id)) %>%
            dplyr::select(hgnc_symbol,
                          everything(),
                          -c(entrezgene_id, ensembl_gene_id)) %>%
            filter(hgnc_symbol!='') %>%
            distinct(hgnc_symbol, .keep_all = TRUE)

        # Store the results in the list
        results_dds.annot[[i]] <- temp_df
        names(results_dds.annot)[i] <- names(results_dds)[i]

    }

    return(results_dds.annot)
}
