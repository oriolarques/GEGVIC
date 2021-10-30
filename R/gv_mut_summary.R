#' @title gv_mut_summary
#'
#' @param muts
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
#' @importFrom maftools read.maf
#' @importFrom maftools plotmafSummary
#' @importFrom maftools oncoplot
#' @import rlang
#' @import dplyr
#'
#' @examples
#' gv_mut_summary(muts = input_gv_module,
#'                metadata = metadata_ge_module,
#'                response = Response,
#'                top_genes = 10,
#'                specific_genes = NULL,
#'                colors = c('black', 'orange'))
#'
#'
gv_mut_summary <- function(muts,
                           metadata,
                           response,
                           top_genes = 10,
                           specific_genes = NULL,
                           colors = c('black' ,'orange')){

    # Process input as MAF file -----------------------------------------------
    maf <- maftools::read.maf(maf = muts,
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

    # Trick to not use quasiquotation with maftools oncoplot to define annotationColor
    ## Get the colname of the response variable so it is quoted
    quoted.resp <- metadata %>%
        dplyr::select(!!response) %>%
        colnames(.)
    ## Get the levels in the response variable
    resp.levels <- metadata %>%
        dplyr::select(!!response) %>%
        dplyr::pull(!!response)
    ## Associate user input colors to an object
    cols <- colors
    # Name the colors vector with the response variable names
    names(cols) <- unique(resp.levels)
    ## Create a list with the named color vector
    cols.list <- list(cols)
    ## Add the quoted named of the response variable as name of the list
    names(cols.list) <- quoted.resp

    # Plot oncoplot
    maftools::oncoplot(maf = maf,
                       top = top_genes,
                       clinicalFeatures = quoted.resp,
                       sampleOrder = samples_order,
                       showTumorSampleBarcodes = TRUE,
                       sortByAnnotation = TRUE,
                       annotationColor = cols.list)
    par(mfrow = c(1,1))



}
