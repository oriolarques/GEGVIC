#' @title ic_plot_comp_samples
#'
#' @description Plots predicted immune cell populations to be compared between
#' sample groups. Among different immune cells, the following are represented:
#' B cells, Macrophages, myeloid Dendritic cells, Neutrophils, Natural Killer,
#' T cell CD4+ and T cell CD8+.
#'
#' @param df Output ot ic_deconv function. A data frame containing the predictions
#' of the six methods included in the immunedeconv package. First column must
#' indicate and be named as cell_type, method and the rest being the numeric
#' prediction for each sample.
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
#' @param points Logical value to decide if points are added to the plot.
#'
#' @return Returns list containing two ggplot objects.
#'
#' @export
#'
#' @import tibble
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @import ggpubr
#' @import rlang
#'
#' @examples
#' tpm <- ic_raw_to_tpm(counts = sample_counts,
#'                      genes_id = 'ensembl_gene_id',
#'                      biomart = ensembl_biomart_GRCh38_p13)
#' ic.pred <- ic_deconv(gene_expression = tpm,
#'                      indications = rep('coad', ncol(tpm)),
#'                      cibersort = NULL,
#'                      tumor = TRUE,
#'                      rmgenes = NULL,
#'                      scale_mrna = TRUE,
#'                      expected_cell_types = NULL)
#' ic_plot_comp_samples(df = ic.pred,
#'                      metadata = sample_metadata,
#'                      response = MSI_status,
#'                      compare = 'wilcox.test',
#'                      p_label = 'p.format',
#'                      colors = c('orange', 'black'),
#'                      points = TRUE)

ic_plot_comp_samples <- function(df,
                                 metadata,
                                 response,
                                 compare = NULL,
                                 p_label = 'p.format',
                                 colors = c('black', 'orange'),
                                 points = TRUE) {

    # Create an object to process the data
    temp_df <- NULL

    # enquote response variable
    response <- rlang::enquo(response)

    # Read cell_type grouping data
    ic_grouping <- GEGVIC:::ic_grouping

    # Create rows with missing data to be added to the predictions data
    ## EPIC and TIMER lack mDendritic, Neutrophils and NK cells
    miss_cat <- as.data.frame(row.names = c('mDendritic_Cell_EPIC',
                                            'Neutrophil_EPIC',
                                            'NK_Cell_TIMER'),
                              matrix(nrow = 3, ncol = (ncol(df) - 2 ))) %>%
        tibble::rownames_to_column('cell_type')

    # Read the results of immune prediction from ic_deconv
    temp_df <- df %>%
        # Reunite cell_type and method columns
        tidyr::unite(data = .,
                     col = 'cell_type',
                     c(cell_type, method),
                     sep = '_') %>%
        # Add rows of missing data created before in miss_cat object
        dplyr::bind_rows(., miss_cat) %>%
        # Join cell_type grouping categories with the ic_grouping data.frame
        dplyr::inner_join(x = .,
                          y = ic_grouping,
                          by = 'cell_type') %>%
        # Transform data into long format
        tidyr::pivot_longer(data = .,
                            cols = c(-cell_type, - grouping),
                            names_to = 'Samples',
                            values_to = 'estimation') %>%
        # Add Samples information with the metadata data.frame
        dplyr::inner_join(x = .,
                          y = metadata,
                          by = 'Samples') %>%
        # Re-separate cell_type and method columns
        tidyr::separate(data = ., col = cell_type, into = c('cell_type', 'method'), sep = '_') %>%
        # Convert columns to factors
        dplyr::mutate(cell_type = as.factor(cell_type),
                      method = as.factor(method),
                      grouping = as.factor(grouping)) %>%
        # Select the necessary columns
        dplyr::select(cell_type, method, grouping, !!response, Samples, estimation) %>%
        # Exclude cell_type classified as Other and CIBERSORT method
        dplyr::filter(grouping != 'Other',
                      method != 'CIBERSORT') %>%
        # Group cell_types by Sample, grouping category and method
        dplyr::group_by(Samples, grouping, method, !!response) %>%
        # Sum the prediction values within different cell_types of the same
        ## category in each patient
        dplyr::summarise(
                         estimation = sum(estimation)) %>%
        # Eliminate duplicated rows
        dplyr::distinct()

    # Plot Between groups comparison per population
    p <- ggplot(temp_df, aes(x = !!response,
                             y = estimation,
                             col = !!response)) +
        # Geometric objects
        geom_violin() +
        geom_boxplot(width = 0.1, outlier.shape = NA) +

        # Define colors
        scale_color_manual(values = colors) +

        # Title
        ggtitle('Immune composition:\nBetween groups comparison per population') +

        # Themes
        theme_bw() +
        theme(
            plot.title = element_text(size = 15, hjust = 0.5, face = 'bold'),
            axis.text.x.bottom = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = 'bottom',
            legend.title = element_text(face='bold', size =12),
            legend.text = element_text(size =12),
            strip.background = element_rect(
                color="black", fill="black", linewidth=1.5, linetype="solid"),
            strip.text = element_text(color = 'white')
        ) +

        # Faceting
        facet_grid(method ~ grouping,
                   scales = 'free_y',
                   switch = 'y')

    # Add p-values for comparisons
    if(is.null(compare) == FALSE){
        p <- p +
            ggpubr::stat_compare_means(method = compare,
                                       label = p_label,
                                       label.y.npc = 0.95,
                                       label.x.npc = 0.3,
                                       show.legend = FALSE)

    }

    # Add points to the plot
    if(points == TRUE){
        p <- p +
            geom_point(alpha = 0.5, position = position_jitter(0.2))
    }

    # Return the plot
    print(p)

}
