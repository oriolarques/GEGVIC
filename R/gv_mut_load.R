#' @title gv_mut_load
#'
#' @description Summarises the total number of mutations per sample and compares
#' different groups of interest.
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
#' @param compare A character string indicating which method to be used for
#' comparing means. Options are 't.test' and 'wilcox.test' for two groups or
#' 'anova' and 'kruskal.test' for more groups. Default value is NULL.
#' @param p_label Character string specifying label type. Allowed values include
#' 'p.signif' (shows the significance levels), 'p.format' (shows the formatted
#' p-value).
#' @param colors Character vector indicating the colors of the different groups
#' to compare. Default values are two: black and orange.
#'
#' @return Returns a ggplot object and a data frame with the data.
#'
#' @export
#'
#' @import dplyr
#' @import rlang
#' @import ggplot2
#' @import ggpubr
#'
#' @examples
#' mut.load <- gv_mut_load(muts = input_gv_module,
#'                         metadata = metadata_ge_module,
#'                         response = Response,
#'                         compare = 'wilcox.test',
#'                         p_label = 'p.format',
#'                         colors = c('black', 'orange'))
#'
gv_mut_load <- function(muts,
                        metadata,
                        response,
                        compare = NULL,
                        p_label = 'p.format',
                        colors = c('black', 'orange')){

    # Get nonsynonymous mutations
    mut.load <- muts %>%
        dplyr::filter(grepl('Missense_Mutation|Nonsense_Mutation', Variant_Classification)) %>%
        # Group by patient
        dplyr::group_by(Tumor_Sample_Barcode) %>%
        # Sum the number of mutations per patient
        dplyr::summarise(mut_load = n()) %>%
        # Join mutations with metadata to get the response groups
        dplyr::right_join(x = .,
                          y = metadata,
                          by = c('Tumor_Sample_Barcode' = 'Samples'))


    # Enquote response variable
    response <-  rlang::enquo(response)

    # Plot mutational load by groups
    p <- ggplot(mut.load, aes(x = !!response,
                              y = mut_load,
                              col = !!response)) +
        # Geometric objects
        geom_violin() +
        geom_boxplot(width = 0.2, outlier.shape = NA) +
        geom_point(alpha = 0.5, position = position_jitter(0.1)) +

        # Define colors
        scale_color_manual(values = colors) +

        # Title
        ggtitle('Mutational Load:\nBetween groups comparison per population') +

        # Themes
        theme_bw() +
        theme(
            plot.title = element_text(size = 15, hjust = 0.5, face = 'bold'),
            axis.text.x.bottom = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = 'bottom'
        )

    # Add p-values for comparisons
    if(is.null(compare) == FALSE){
        p <- p +
            ggpubr::stat_compare_means(method = compare,
                                       label = p_label,
                                       label.y.npc = 0.95,
                                       label.x.npc = 0.3,
                                       show.legend = FALSE)

    }

    # Print the plot
    print(p)
    # Return the table
    return(mut.load)

}


