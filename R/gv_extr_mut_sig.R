#' @title gv_extr_mut_sig
#'
#' @description
#'
#' @param results
#' @param ids_samples
#'
#' @return
#'
#' @export
#'
#' @import tibble
#' @import tidyr
#' @import dplyr
#'
#' @examples
#' mut.sbs <- input_gv_module %>% dplyr::filter(Variant_type == 'SNP')
#' sigs.sbs.input <- deconstructSigs::mut.to.sigs.input(
#'                                           mut.ref = mut.sbs,
#'                                           sample.id = 'Tumor_Sample_Barcode',
#'                                           chr = 'Chromosome',
#'                                           pos = 'Start_Position',
#'                                           ref = 'Reference_Allele',
#'                                           alt = 'Tumor_Seq_Allele2',
#'                                           bsg = BSgenome.Hsapiens.UCSC.hg19,
#'                                           sig.type = 'SBS')
#' ids_samples <- unique(input_gv_module$Tumor_Sample_Barcode)
#' results_sbs <- sapply(ids_samples,
#'                       function(x) {
#'                          deconstructSigs::whichSignatures(
#'                          tumor.ref = sigs.sbs.input,
#'                          signatures.ref = signatures.cosmic,
#'                                  sample.id = x,
#'                                  contexts.needed = TRUE,
#'                                  tri.counts.method = 'default')#'
#'
#' results_sbs.extr <- gv_extr_mut_sig(results_sbs)
#'
gv_extr_mut_sig <- function(results,
                            ids_samples){
    # format results output ---------------------------------------------------
    results.sig <- as.data.frame(results)
    results.sig <- results.sig[1,]  # filter signature weigth results
    results.sig <- as.data.frame(unlist(results.sig))
    names(results.sig) <- 'results'

    # Convert rownames to column Sample ---------------------------------------
    results.sig <- results.sig %>%
        tibble::rownames_to_column('samples')

    names(results.sig) <- c('samples', 'Value')

    # Separate the sample colum into sample and singature ---------------------
    results.sig <- tidyr::separate(data = results.sig,
                                   col = samples,
                                   into = c("Samples", "Signature"),
                                   sep = ".weights.")

    # Rename samples with the original sample name
    results.sig$Samples <- as.factor(results.sig$Samples)
    levels(results.sig$Samples) <- ids_samples


    # Fix the numeration of the signature below 10 ----------------------------
    results.sig <- results.sig %>%
        # change signature1 for singature01
        dplyr::mutate(Signature = gsub(pattern = '^([aA-zZ]+\\.*)([0-9])$',
                                       replacement = '\\10\\2',
                                       x = Signature))

    return(results.sig)

}
