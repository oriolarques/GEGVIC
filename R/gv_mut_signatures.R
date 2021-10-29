#' Title
#'
#' @param muts
#' @param metadata
#' @param response
#' @param gbuild
#' @param mut_sigs
#' @param tri.counts.method
#' @param colors
#'
#' @return
#'
#' @export
#'
#' @import dplyr
#' @import deconstructSigs
#' @import ggplot2
#' @import ggpubr
#' @import tidyr
#' @import tibble
#' @import pheatmap
#' @import BSgenome.Hsapiens.UCSC.hg19
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @import ggplotify
#'
#' @examples
#'
#' gv_mut_signatures(muts = input_gv_module,
#'                   metadata = metadata_ge_module,
#'                   response = Response,
#'                   gbuild = BSgenome.Hsapiens.UCSC.hg19,
#'                   mut_sigs = signatures.cosmic,
#'                   tri.counts.method = 'default',
#'                   colors = c('black', 'orange'))
#'
gv_mut_signatures <- function(muts,
                              metadata,
                              response,
                              gbuild,
                              mut_sigs = signatures.cosmic,
                              tri.counts.method = 'default',
                              colors = c('black', 'orange')) {

    # Enquote response variable
    response <- rlang::enquo(response)

    # Split the mutations input into SNP and DNP ------------------------------
    mut.sbs <- muts %>%
        dplyr::filter(Variant_type == 'SNP')

    mut.dbs <- muts %>%
        dplyr::filter(Variant_type == 'DNP')

    # Create deconstructSigs inputs -------------------------------------------
    # If single base substitutions are provided
    if (nrow(mut.sbs != 0)) {
        sigs.sbs.input <- deconstructSigs::mut.to.sigs.input(
            mut.ref = mut.sbs,
            sample.id = 'Tumor_Sample_Barcode',
            chr = 'Chromosome',
            pos = 'Start_Position',
            ref = 'Reference_Allele',
            alt = 'Tumor_Seq_Allele2',
            bsg = gbuild,
            sig.type = 'SBS')
    }
    # If dobulet base substitutions are provided
    if (nrow(mut.dbs != 0)) {
        sigs.dbs.input <- deconstructSigs::mut.to.sigs.input(
            mut.ref = mut.dbs,
            sample.id = 'Tumor_Sample_Barcode',
            chr = 'Chromosome',
            pos = 'Start_Position',
            ref = 'Reference_Allele',
            alt = 'Tumor_Seq_Allele2',
            bsg = gbuild,
            sig.type = 'DBS')
    }

    # generate ids for all samples ------------------------------------------------
    ids_samples <- unique(muts$Tumor_Sample_Barcode)

    # get mutational signature predictions for all samples ------------------------
    # If single base substitutions are provided
    if (nrow(mut.sbs != 0)) {
        results_sbs <- sapply(ids_samples,
                              function(x) {
                                  deconstructSigs::whichSignatures(
                                      tumor.ref = sigs.sbs.input,
                                      signatures.ref = mut_sigs,
                                      sample.id = x,
                                      contexts.needed = TRUE,
                                      tri.counts.method = tri.counts.method)
                              })
    }
    # If dobulet base substitutions are provided
    if (nrow(mut.dbs != 0)) {
        results_dbs <- sapply(ids_samples,
                              function(x) {
                                  deconstructSigs::whichSignatures(
                                      tumor.ref = sigs.dbs.input,
                                      signatures.ref = mut_sigs,
                                      sample.id = x,
                                      contexts.needed = TRUE,
                                      tri.counts.method = tri.counts.method)
                              })
    }

    # Analyze results --------------------------------------------------------
    # If single base substitutions are provided
    if (nrow(mut.sbs != 0)) {
        # Extract results from whichSignatures function
        results_sbs.extr <- gv_extr_mut_sig(results = results_sbs,
                                            ids_samples = ids_samples) %>%
            # Join predicted mutational signature results with metadata
            dplyr::left_join(x = .,
                             y = metadata_ge_module,
                             by = c('Samples')) %>%
            # Round predicted mutational signature contribution
            dplyr::mutate(Value = round(x = Value, digits = 2))
    }
    # If dobulet base substitutions are provided
    if (nrow(mut.dbs != 0)) {
        # Extract results from whichSignatures function
        results_dbs.extr <- gv_extr_mut_sig(results = results_dbs,
                                            ids_samples = ids_samples) %>%
            # Join predicted mutational signature results with metadata
            dplyr::left_join(x = .,
                             y = metadata_ge_module,
                             by = c('Samples')) %>%
            # Round predicted mutational signature contribution
            dplyr::mutate(Value = round(x = Value, digits = 2))
    }



    # Plot results ------------------------------------------------------------
    ## Barplot  --------------------------------------------------------------
    bar.plot <- ggplot(results_sbs.extr, aes(x = Samples,
                                             y = Value,
                                             fill = as.factor(Signature))) +

        # Geometric objects
        geom_bar(stat = 'identity') +

        # Define fill colors using the Set1 palette from ggpubr package
        scale_fill_manual(values = ggpubr::get_palette(palette = 'Set1',
                                                       k = length(unique(
                                                           results_sbs.extr$Value
                                                       )))) +
        # Expand columns to fill margins
        scale_y_continuous(expand = c(0,0)) +

        # Title and labs
        ggtitle('Mutational signature predictions') +
        labs(fill = 'Signatures') +

        # Themes
        theme_linedraw() +
        theme(
            plot.title = element_text(size = 15, hjust = 0.5, face = 'bold'),
            #axis.text.x.bottom = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = 5, angle = 45, hjust = 1)
        ) +

        # Faceting
        facet_wrap(facets = vars(!!response),
                   scales = 'free_x')


    ## Heatmap  ------------------------------------------------------------
    # Format signature predictions object in a wide format: Pivot wider
    wide.results_sbs.extr <- results_sbs.extr %>%
        dplyr::select(Samples, Signature, Value) %>%
        tidyr::pivot_wider(id_cols = Signature,
                           names_from = Samples,
                           values_from = Value) %>%
        tibble::column_to_rownames('Signature')



    # Format the response variable from metadata
    pheat.meta <- metadata %>%
        dplyr::select(Samples, !!response) %>%
        dplyr::arrange(!!response) %>%
        tibble::column_to_rownames('Samples')

    # Define response level group colors in a list
    temp_color <- colors

    resp.levels <- metadata %>%
        dplyr::select(!!response) %>%
        dplyr::pull(!!response)

    names(temp_color) <- unique(resp.levels)

    # Get the quoted name of the response variable
    quoted.resp <- metadata %>%
        dplyr::select(!!response) %>%
        colnames(.)

    pheat.anno.color <- list(temp_color)
    # Name the list
    names(pheat.anno.color) <- quoted.resp

    # Plot pheatmap
    heat.map <- pheatmap(as.matrix(wide.results_sbs.extr[,
                                                         order(match(colnames(wide.results_sbs.extr),
                                                                     rownames(pheat.meta)))]),
                         color = ggpubr::get_palette(palette = 'Purples', k = 10),
                         scale = 'none',
                         cluster_rows = FALSE,
                         cluster_cols = FALSE,
                         annotation_col = pheat.meta,
                         annotation_colors = pheat.anno.color,
                         silent = TRUE)

    # Merge the resulting plots -----------------------------------------------
    mut.plot.list <- list(mut_sig_barplot = bar.plot,
                          mut_sig_heatmap = ggplotify::as.ggplot(heat.map))

    return(mut.plot.list)

}
