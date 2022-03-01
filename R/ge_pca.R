#' @title ge_pca
#'
#' @description Creates a PCA plot from gene expression data using the VST
#' transformation from the plotPCA function from the DESeq2 package.
#'
#' @param counts Data frame that contains gene expression data as raw counts.
#' @param genes_id Name of the column that contains gene identifiers. Should be
#' one of the following:'entrezgene_id', 'ensembl_gene_id' or 'hgnc_symbol'.
#' @param metadata Data frame that contains supporting variables to the data.
#' @param response Unquoted name of the variable indicating the groups to analyse.
#' @param design Variables in the design formula in the form of: 'Var1 + Var2 + ... Var_n'.
#' @param colors Character vector indicating the colors of the different groups
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
#' @import rlang
#'
#' @examples
#' ge_pca(counts = input_ge_module,
#'        genes_id = 'entrezgene_id',
#'        metadata = metadata_ge_module,
#'        response = Response,
#'        design = 'Response',
#'        colors = c('black', 'orange'))
#'
ge_pca <- function(counts,
                   genes_id,
                   metadata,
                   response,
                   design,
                   colors = c('black', 'orange')) {

    # Preprocess counts data
    counts <- preprocess_ge_counts(counts = counts,
                                   genes_id = genes_id)

    # Preprocess metadata
    metadata <- preprocess_ge_meta(metadata = metadata)

    # Create DESeq2Dataset object
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
                                          colData = metadata,
                                          design = formula(paste('~', design,
                                                                 collapse = " ")))
    # Generate  normalized counts
    dds <- DESeq2::estimateSizeFactors(dds)

    # Transform normalized counts for data visualization
    vsd <- DESeq2::vst(dds, blind=FALSE)

    ## Enquote response variable
    response <- rlang::enquo(response)
    # Trick to not use quasiquotation in plotPCA intgroup argument
    ## Get the colname of the response variable so it is quoted
    quoted.resp <- metadata %>%
        dplyr::select(!!response) %>%
        colnames(.)

    # plot PCA
    DESeq2::plotPCA(vsd, intgroup = quoted.resp) +
                scale_color_manual(values = colors) +
                labs(title = 'Principal Component Analysis') +
                theme_bw() +
                theme(plot.title = element_text(hjust = 0.5))

}
