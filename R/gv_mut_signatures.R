#' @title gv_mut_signatures
#'
#' @description Predicts the contribution of known mutational processes to the samples.
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
#' - Variant_Type: Type of mutation. Can be: 'SNP' (Single nucleotide polymorphism),
#' 'DNP' (Double nucleotide polymorphism), 'INS' (Insertion), 'DEL' (Deletion).
#' - Tumor_Sample_Barcode: Sample name.
#' @param metadata Data frame that contains supporting variables to the data.
#' @param response Unquoted name of the variable indicating the groups to analyse.
#' @param gbuild Version of the genome to work with. It can be one of the following:
#' - ‘BSgenome.Hsapiens.UCSC.hg19’
#' - ‘BSgenome.Hsapiens.UCSC.hg38’
#' - ‘BSgenome.Mmusculus.UCSC.mm10’
#' - ‘BSgenome.Mmusculus.UCSC.mm39’
#' @param mut_sigs Mutational signature matrices containing the frequencies of
#' all nucleotide changes per signature need to be indicated. GEGVIC contains the
#' matrices from COSMIC for single and double base substitutions. To choose one,
#' the user has to indicate ’COSMIC_v{XX}_{YY}BS_GRCh{ZZ}’ in the mut_sigs argument.
#' The XX is the version, that can be v2 or v3.2. YY indicates if mutations are
#' single (S) or double (D) base substitutions, while the ZZ is for the genome
#' assembly, either GRCh37 or GRCh38 for human data and mm9 or mm10 for mouse data.
#' @param tri.counts.method Normalization method. Needs to be set to either:
#' - 'default' – no further normalization.
#' - 'exome' – normalized by number of times each trinucleotide context is
#' observed in the exome.
#' - 'genome' – normalized by number of times each trinucleotide context is
#' observed in the genome.
#' - 'exome2genome' – multiplied by a ratio of that trinucleotide's occurence in
#' the genome to the trinucleotide's occurence in the exome.
#' - 'genome2exome' – multiplied by a ratio of that trinucleotide's occurence in
#' the exome to the trinucleotide's occurence in the genome.
#' - data frame containing user defined scaling factor – count data for each
#' trinucleotide context is multiplied by the corresponding value given in the
#' data frame.
#' @param colors Character vector indicating the colors of the different groups
#' to compare. Default values are two: black and orange.
#' @param col.names Logical value to determine if tumour are shown in plots.
#'
#' @return Returns ggplot objects and the results table in a form of a data frame.
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
#' mut.sigs <- gv_mut_signatures(muts = sample_mutations,
#'                               metadata = sample_metadata,
#'                               response = MSI_status,
#'                               gbuild = 'BSgenome.Hsapiens.UCSC.hg38',
#'                               mut_sigs = 'COSMIC_v2_SBS_GRCh38',
#'                               tri.counts.method = 'default',
#'                               colors = c('orange', 'black'),
#'                               col.names = TRUE)
#'
gv_mut_signatures <- function(muts,
                              metadata,
                              response,
                              gbuild = 'BSgenome.Hsapiens.UCSC.hg19',
                              mut_sigs = 'COSMIC_v2_SBS_GRCh37',
                              tri.counts.method = 'default',
                              colors = c('black', 'orange'),
                              col.names = TRUE) {

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
            dplyr::filter(Variant_Type == 'SNP')
        # Define sig.type as SBS
        sig_type <- 'SBS'

    } else if (grepl('DBS', eval.mut.input, ignore.case = TRUE) == TRUE) {
        # Filter mutations of DNP type
        mut.filt <- muts %>%
            dplyr::filter(Variant_Type == 'DNP')
        # Define sig.type as DBS
        sig_type <- 'DBS'

    } else {
        # Filter mutations of INS or DEL type
        mut.filt <- muts %>%
            dplyr::filter(Variant_Type %in% c('INS', 'DEL'))
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

    # Analyze results ---------------------------------------------------------
    # Extract results from whichSignatures function
    results.extr <- GEGVIC::gv_extr_mut_sig(results = results,
                                            ids_samples = ids_samples) %>%
        # Join predicted mutational signature results with metadata
        dplyr::left_join(x = .,
                         y = metadata,
                         by = c('Samples')) %>%
        # Round predicted mutational signature contribution
        dplyr::mutate(Value = round(x = Value, digits = 2))

    # Filter top 4 signatures for barplot
    top.results.extr <- results.extr %>%
        dplyr::group_by(Samples) %>%
        dplyr::top_n(n = 4, wt = Value) %>%
        droplevels()



    # Plot results ------------------------------------------------------------
    ## Barplot  ---------------------------------------------------------------
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
        theme_bw() +
        theme(
            plot.title = element_text(size = 15, hjust = 0.5, face = 'bold'),
            #axis.text.x.bottom = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = 8, angle = 45, hjust = 1, face = 'bold'),
            strip.background = element_rect(
                color="black", fill="black", size=1.5, linetype="solid"),
            strip.text = element_text(color = 'white')
        ) +

        # Faceting
        facet_wrap(facets = vars(!!response),
                   scales = 'free_x')

    ## Eliminate sample names if the user decides so
    if(col.names == FALSE){
        bar.plot <- bar.plot + theme(axis.text.x = element_blank())
    }

    ## Heatmap  ---------------------------------------------------------------
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
                         show_colnames = col.names,
                         scale = 'none',
                         cluster_rows = FALSE,
                         cluster_cols = FALSE,
                         annotation_col = pheat.meta,
                         annotation_colors = pheat.anno.color,
                         main = 'Mutational signature predictions per sample',
                         silent = TRUE)

    # Merge the resulting plots -----------------------------------------------
    print(bar.plot)
    print(ggplotify::as.ggplot(heat.map))

    return(results.extr)

}
