#' @title ge_pca
#'
#' @description
#'
#' @param counts Data frame that contains gene expression data as raw counts.
#' @param genes_id Name of the column that contains gene identifiers. Should be
#' one of the following:'entrez_gene_id', 'ensembl_gene_id' or 'hgnc_symbol'.
#' @param metadata Data frame that contains supporting variables to the data.
#' @param design Variables in the design formula in the form of: 'Var1 + Var2 + ... Var_n'.
#' @param colors  Character vector indicating the colors of the different groups
#' to compare. Default values are two: black and orange.
#'
#' @return Returns a ggplot object.
#'
#' @export
#'
#' @import DESeq2
#' @import dplyr
#' @import ggplot2
#' @import tibble
#'
#' @examples
#'
ge_pca <- function(counts,
                   genes_id,
                   metadata,
                   design,
                   colors = c('black', 'orange')) {

    # Preprocess counts data
    counts <- preprocess_ge_counts(counts = counts,
                                   genes_id = genes_id)

    # Get patient ID's as rownames
    metadata <- preprocess_ge_meta(metadata = metadata)

    # Create DESeq2Dataset object
    dds <- DESeqDataSetFromMatrix(countData = counts,
                                  colData = metadata,
                                  design = formula(paste('~', design, collapse = " ")))
    # Generate  normalized counts
    dds <- estimateSizeFactors(dds)

    # Transform normalized counts for data visualization
    vsd <- vst(dds, blind=FALSE)

    # plot PCA
    plotPCA(vsd, intgroup="Response") +
        scale_color_manual(values = colors) +
        labs(title = 'Principal Component Analysis') +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))

}
