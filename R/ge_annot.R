#' Title
#'
#' @param results_dds
#' @param genes_id Name of the column that contains gene identifiers. Should be
#' one of the following:'entrez_gene_id', 'ensembl_gene_id' or 'hgnc_symbol'.
#'
#' @return
#' @export
#'
#' @examples

ge_annot <- function(results_dds,
                     genes_id){

    # Create a data.frame to store temporarilly data
    temp_df <- NULL

    genes_id <- enquo(genes_id)

    # Iterate over the differential gene expression results list
    for(i in seq_along(results_dds)){

        # Get the results as a data frame
        temp_df <- as.data.frame(results_dds[i]) %>%
            # Rownames to column with the name indicated in the genes_id parameter
            rownames_to_column(!!genes_id) %>%
            mutate(!!genes_id := as.character(!!genes_id)) %>%
            inner_join(x = .,
                       y = ensembl_biomart_GRCh38_p13 %>%
                           mutate(entrezgene_id = as.character(entrezgene_id)),
                       by = gene_id) %>%
            dplyr::select(hgnc_symbol,
                          everything(),
                          -c(entrezgene_id, ensembl_gene_id)) %>%
            filter(hgnc_symbol!='') %>%
            distinct(hgnc_symbol, .keep_all = TRUE) %>%
            column_to_rownames('hgnc_symbol')

        #


    }

}
