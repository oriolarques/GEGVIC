#' @title gv_mut_summary
#'
#' @param df
#' @param metadata
#' @param response
#' @param top_genes
#' @param specific_genes
#' @param colors
#'
#' @return
#'
#' @export
#'
#' @import maftools
#' @import rlang
#' @import dplyr
#'
#' @examples
gv_mut_summary <- function(df,
                           metadata,
                           response,
                           top_genes = 10,
                           specific_genes = NULL,
                           colors = c('black' ,'orange')){

    # Process input as MAF file -----------------------------------------------
    maf <- maftools::read.maf(maf = df,
                              clinicalData = metadata %>%
                                  dplyr::rename('Tumor_Sample_Barcode' = 'Samples'))

    # Plot MAF Summary --------------------------------------------------------
    maftools::plotmafSummary(maf = maf,
                             addStat = 'median',
                             titvRaw = FALSE,
                             top = top_genes)
    par(mfrow = c(1,1))

    # Oncoplot ----------------------------------------------------------------

    # Enquote response variable
    response <-  rlang::enquo(response)



    # Determine the order of the samples to separate between groups in the plot
    samples_order <- metadata %>%
        dplyr::arrange(!!response) %>%
        dplyr::pull(Samples)


    x <- metadata %>%
        dplyr::select(!!response) %>%
        colnames(.)

    y <- metadata %>%
        dplyr::select(!!response) %>%
        dplyr::pull(!!response)

    cols <- colors
    names(cols) <- unique(y)

    z <- list(cols)
    names(z) <- x

    # Plot oncoplot
    maftools::oncoplot(maf = maf,
                       top = top_genes,
                       clinicalFeatures = x,
                       sampleOrder = samples_order,
                       showTumorSampleBarcodes = TRUE,
                       sortByAnnotation = TRUE,
                       annotationColor = z)
    par(mfrow = c(1,1))



}
