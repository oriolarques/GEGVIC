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
#' counts <- preprocess_ge_counts(counts = sample_counts,
#'                                genes_id = 'ensembl_gene_id')
#'
preprocess_ge_counts <- function(counts,
                                 genes_id) {

    # Make sure the input does not has rownames
    counts <- counts
    rownames(counts) <- c()

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
#' @param counts Data frame that contains gene expression data as raw counts.
#'
#' @return Generates a data frame.
#'
#' @export
#'
#' @import dplyr
#' @import tibble
#'
#' @examples
#' metadata <- preprocess_ge_meta(metadata = sample_metadata,
#'                                counts = sample_counts[,-1])

preprocess_ge_meta <- function(metadata,
                               counts = NULL) {

    # Make sure the input does not has rownames
    metadata <- metadata
    rownames(metadata) <- c()

    # Reorder the metadata so the Samples are in the same order as in counts
    if(is.null(counts) == FALSE){
        metadata <- metadata %>%
            dplyr::arrange(match(Samples, colnames(counts))) %>%
            dplyr::mutate_all(., as.factor) %>%
            # Get patient ID's as rownames
            tibble::column_to_rownames('Samples')
    } else {
        metadata <- metadata <- metadata %>%
            dplyr::mutate_all(., as.factor) %>%
            # Get patient ID's as rownames
            tibble::column_to_rownames('Samples')

    }

    return(metadata)

}
