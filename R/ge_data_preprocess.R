#' @title preprocess_ge_counts
#'
#' @description Prepares RNA-seq raw counts for posterior analysis using DESeq2.
#'
#'
#' @param counts Data frame that contains gene expression data as raw counts.
#' @param genes_id Name of the column that contains gene identifiers. Should be
#' one of the following:'entrezgene_id', 'ensembl_gene_id' or 'hgnc_symbol'.
#'
#' @return Generates a data frame.
#'
#' @export
#'
#' @import dplyr
#' @import tibble
#'
#' @examples
#' counts <- preprocess_ge_counts(counts = counts, genes_id = genes_id)
#'
preprocess_ge_counts <- function(counts,
                                 genes_id) {
    # Get gene symbols as rownames
    counts <- counts %>%
        tibble::column_to_rownames(genes_id)

    return(counts)

}


#' @title preprocess_ge_meta
#'
#' @description Prepares metadata file for posterior analysis using DESeq2.
#'
#' @param metadata Data frame that contains supporting variables to the data.
#'
#' @return Generates a data frame.
#'
#' @export
#'
#' @import dplyr
#' @import tibble
#'
#' @examples
#' metadata <- preprocess_ge_meta(metadata = metadata)

preprocess_ge_meta <- function(metadata) {
    # Get patient ID's as rownames
    metadata <- metadata %>%
        dplyr::mutate_all(., as.factor) %>%
        tibble::column_to_rownames('Samples')

    return(metadata)

}
