#' @title gv_mut_summary
#'
#' @description Given a table that contains genetic variants from samples of
#' interest, it summarises different aspects of such mutations and draws an
#' oncoplot.
#'
#' @param muts Data frame containing genetic variations. Necessary columns must
#' have the following names:
#' - Hugo_Symbol: Gene symbol from HGNC.
#' - Chromosome: Affected chromosome.
#' - Start_Position: Mutation start coordinate.
#' - End_Position: Mutation end coordinate.
#' - Reference_Allele: The plus strand reference allele at this position.
#' Includes the deleted sequence for a deletion or "-" for an insertion.
#' - Tumor_Seq_Allele2: Tumor sequencing discovery allele.
#' - Variant_Classification: Translational effect of variant allele. Can be one
#' of the following: Frame_Shift_Del, Frame_Shift_Ins, In_Frame_Del,
#' In_Frame_Ins, Missense_Mutation, Nonsense_Mutation, Silent, Splice_Site,
#' Translation_Start_Site, Nonstop_Mutation, RNA, Targeted_Region.
#' - Variant_type: Type of mutation. Can be: 'SNP' (Single nucleotide polymorphism),
#' 'DNP' (Double nucleotide polymorphism), 'INS' (Insertion), 'DEL' (Deletion).
#' - Tumor_Sample_Barcode: Sample name.
#' @param metadata Data frame that contains supporting variables to the data.
#' @param response Unquoted name of the variable indicating the groups to analyse.
#' @param top_genes Number of genes to be analysed in the mutational summary.
#' @param specific_genes Genes that will be plotted in the oncoplot.
#' @param colors Character vector indicating the colors of the different groups
#' to compare. Default values are two: black and orange.
#'
#' @return Prints two plots: A summary of samples mutations and an oncoplot.
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
                       genes = specific_genes,
                       sampleOrder = samples_order,
                       showTumorSampleBarcodes = TRUE,
                       sortByAnnotation = TRUE,
                       annotationColor = cols.list)
    par(mfrow = c(1,1))



}
