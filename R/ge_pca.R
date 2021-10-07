#' Title
#'
#' @param data Data frame that contains the gene expression data
#' @param gs
#' @param metadata Data frame that contains supporting variables to the data
#' such as sample groups, treatments, etc.
#' @param group Unquoted name of the column in the metadata where the grouping factor is
#' located.
#' @param colors Character vector indicating the colors of the different groups
#' to compare. Default values are black and orange

#'
#' @return
#' @export
#'
#' @import DESeq2
#' @import dplyr
#' @import ggplot2
#'
#' @examples
#'
ge_pca <- function(data,
                   gs,
                   metadata,
                   group,
                   colors = c('black', 'orange')) {

    group <- enquo(group)

    # Get gene symbols as rownames
    data <- data %>%
        column_to_rownames(gs)

    # Get patient ID's as rownames
    metadata <- metadata %>%
        mutate(!!group := as.factor(!!group)) %>%
        column_to_rownames('Samples')

    print(metadata)
    # Create DESeq2Dataset object
    dds <- DESeqDataSetFromMatrix(countData = data,
                                  colData = metadata,
                                  design = formula(paste('~', group, collapse = " ")))
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
