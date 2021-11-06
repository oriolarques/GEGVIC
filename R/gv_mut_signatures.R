#' Title
#'
#' @param muts
#' @param metadata
#' @param response
#' @param gbuild
#' @param mut_sigs can be COSMIC_v2_SBS_GRCh37, COSMIC_v2_SBS_GRCh38 ...
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
#' @import ggplotify
#' @import rlang
#'
#' @examples
#'
#' gv_mut_signatures(muts = input_gv_module,
#'                   metadata = metadata_ge_module,
#'                   response = Response,
#'                   gbuild = 'BSgenome.Hsapiens.UCSC.hg19',
#'                   mut_sigs = 'COSMIC_v2_SBS_GRCh37',
#'                   tri.counts.method = 'default',
#'                   colors = c('black', 'orange'))
#'
gv_mut_signatures <- function(muts,
                              metadata,
                              response,
                              gbuild = 'BSgenome.Hsapiens.UCSC.hg19',
                              mut_sigs = COSMIC_v2_SBS_GRCh37,
                              tri.counts.method = 'default',
                              colors = c('black', 'orange')) {

    # Enquote response variable
    response <- rlang::enquo(response)

    # Load genomic build
    library(gbuild, character.only = TRUE)

    # Check the type of mutations to use --------------------------------------

    ## Get the name of the mutational signatures files chosen by the user
    eval.mut.input <- substitute(mut_sigs)

    # If the name of the mutational singatures contains SBS
    if (grepl('SBS', eval.mut.input, ignore.case = TRUE) == TRUE) {
        # Filter mutations of SNP type
        mut.filt <- muts %>%
            dplyr::filter(Variant_type == 'SNP')
        # Define sig.type as SBS
        sig_type <- 'SBS'

    } else if (grepl('DBS', eval.mut.input, ignore.case = TRUE) == TRUE) {
        # Filter mutations of DNP type
        mut.filt <- muts %>%
            dplyr::filter(Variant_type == 'DNP')
        # Define sig.type as DBS
        sig_type <- 'DBS'

    } else {
        # Filter mutations of INS or DEL type
        mut.filt <- muts %>%
            dplyr::filter(Variant_type %in% c('INS', 'DEL'))
        # Define sig.type as SBS
        sig_type <- 'ID'

    }


    # Create deconstructSigs inputs -------------------------------------------
    sigs.input <- deconstructSigs::mut.to.sigs.input(
        mut.ref = mut.filt,
        sample.id = 'Tumor_Sample_Barcode',
        chr = 'Chromosome',
        pos = 'Start_Position',
        ref = 'Reference_Allele',
        alt = 'Tumor_Seq_Allele2',
        bsg = get(noquote(gbuild)),
        sig.type = sig_type)


    # generate ids for all samples --------------------------------------------
    ids_samples <- unique(muts$Tumor_Sample_Barcode)

    # get mutational signature predictions for all samples --------------------
    results <- sapply(ids_samples,
                      function(x) {
                          deconstructSigs::whichSignatures(
                              tumor.ref = sigs.input,
                              signatures.ref = as.data.frame(get(noquote(mut_sigs))),
                              sample.id = x,
                              contexts.needed = TRUE,
                              tri.counts.method = tri.counts.method)
                      })

    # Analyze results --------------------------------------------------------
    # Extract results from whichSignatures function
    results.extr <- GEGVIC::gv_extr_mut_sig(results = results,
                                            ids_samples = ids_samples) %>%
        # Join predicted mutational signature results with metadata
        dplyr::left_join(x = .,
                         y = metadata_ge_module,
                         by = c('Samples')) %>%
        # Round predicted mutational signature contribution
        dplyr::mutate(Value = round(x = Value, digits = 2))

    # Filter top 4 signatures for barplot
    top.results.extr <- results.extr %>%
        dplyr::group_by(Samples) %>%
        dplyr::top_n(n = 4, wt = Value) %>%
        droplevels()



    # Plot results ------------------------------------------------------------
    ## Barplot  --------------------------------------------------------------
    bar.plot <- ggplot(top.results.extr, aes(x = Samples,
                                             y = Value,
                                             fill = as.factor(Signature))) +

        # Geometric objects
        geom_bar(stat = 'identity') +

        # Define fill colors using the Set1 palette from ggpubr package
        scale_fill_manual(values = ggpubr::get_palette(palette = 'simpsons',
                                                       k = length(unique(
                                                           top.results.extr$Signature
                                                       )))) +
        # Expand columns to fill margins
        scale_y_continuous(expand = c(0,0)) +

        # Title and labs
        ggtitle('Top 4 Mutational signature predictions per sample') +
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
    wide.results.extr <- results.extr %>%
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
    heat.map <- pheatmap(as.matrix(wide.results.extr[,
                                                         order(match(colnames(wide.results.extr),
                                                                     rownames(pheat.meta)))]),
                         color = ggpubr::get_palette(palette = 'Purples', k = 10),
                         scale = 'none',
                         cluster_rows = FALSE,
                         cluster_cols = FALSE,
                         annotation_col = pheat.meta,
                         annotation_colors = pheat.anno.color,
                         main = 'Mutational signature predictions per sample',
                         silent = TRUE)

    # Merge the resulting plots -----------------------------------------------
    mut.plot.list <- list(mut_sig_barplot = bar.plot,
                          mut_sig_heatmap = ggplotify::as.ggplot(heat.map))

    return(mut.plot.list)

}
