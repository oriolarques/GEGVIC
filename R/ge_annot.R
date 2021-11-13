#' @title ge_annot
#'
#' @description Associates Entrez_gene_id (NCBI-ID) or ENSEMBL_gene_id to
#' HGNC symbols (Hugo Gene Nomenclature Committee nomenclature).
#'
#' @param results_dds The output of the ge_diff_exp function. List containing
#' data frames of differential gene expression results between the different groups.
#' @param genes_id Name of the column that contains gene identifiers. Should be
#' one of the following: 'ensembl_gene_id', 'entrezgene_id' or 'hgnc_symbol'.
#' @param biomart Data frame containing a biomaRt query with the following
#' attributes: ensembl_gene_id, hgnc_symbol, entrezgene_id, transcript_length,
#' refseq_mrna. In the case of mus musculus data, external_gene_name must be
#' obtained and then change the column name for hgnc_symbol. Uploaded biomaRt
#' queries in GEGVIC: 'ensembl_biomartGRCh37', ensembl_biomartGRCh38_p13' and
#' 'ensembl_biomartGRCm38_p6', 'ensembl_biomartGRCm39'.
#'
#' @return Returns a list of data frames containing differential gene expression
#' results, one for each level comparison where genes have been annotated using
#' HGNC symbols.
#'
#' @export
#'
#' @import DESeq2
#' @import dplyr
#' @import tibble
#'
#' @examples
#' results_dds <- ge_diff_exp(counts = input_ge_module,
#'                            genes_id = 'entrezgene_id',
#'                            metadata = metadata_ge_module,
#'                            design = 'Response',
#'                            ref_level = c('Response', 'Non_Responders'),
#'                            shrink = 'apeglm')
#' annot.res <- ge_annot(results_dds = results_dds,
#'                       genes_id = 'entrezgene_id',
#'                       biomart = ensembl_biomart_GRCh38_p13)

ge_annot <- function(results_dds,
                     genes_id,
                     biomart){

    # Create a data.frame to store temporarilly data
    temp_df <- NULL

    # Create a list to store the annotated results
    results_dds.annot <- list()

    # Iterate over the differential gene expression results list
    for (i in seq_along(results_dds)) {

        # Get the results as a data frame
        # First evaluate if genes are annotated already as hgnc_symbol
        if (genes_id == 'hgnc_symbol'){
            # If so:
            temp_df <- as.data.frame(results_dds[[i]]) %>%
                # Rownames to column with the name indicated in the genes_id parameter
                tibble::rownames_to_column(genes_id) %>%
                # Filter those missing gene symbols
                dplyr::filter(hgnc_symbol != '') %>%
                # Remove duplicated genes
                dplyr::distinct(hgnc_symbol, .keep_all = TRUE)

        } else {
            # If genes are identified as entrezgene_id or ensembl_gene_id:
            temp_df <- as.data.frame(results_dds[[i]]) %>%
                # Rownames to column with the name indicated in the genes_id parameter
                tibble::rownames_to_column(genes_id) %>%
                # Join the data frame with the GRCh38_p13 biomaRt table stored in data
                dplyr::inner_join(x = .,
                           y = biomart %>%
                               mutate(entrezgene_id = as.character(entrezgene_id)),
                           by = genes_id) %>%
                # From all annotation columns keep only hgnc symbol column
                dplyr::select(hgnc_symbol,
                              everything(),
                              -c(entrezgene_id, ensembl_gene_id,
                                 transcript_length, refseq_mrna)) %>%
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
