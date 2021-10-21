#' @title ic_plot_comp_samples
#'
#' @description
#'
#' @param df
#' @param metadata
#' @param response
#' @param compare
#' @param p_label
#' @param colors
#'
#' @return Returns a ggplot object.
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
#'
#'
#'
#'
#'
ic_plot_comp_samples <- function(df,
                                 metadata,
                                 response, # unquoted name of grouping variable
                                 compare = NULL, # compare method ggpubr: t.test, wilcox.test, anova, kruskal.test
                                 p_label = 'p.format',
                                 colors = c('black', 'orange')) {

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
                          y = metadata_ge_module,
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
        geom_point(alpha = 0.5, position = position_jitter(0.2)) +

        # Define colors
        scale_color_manual(values = colors) +

        # Title
        ggtitle('Immune composition:\nBetween groups comparison per population') +

        # Themes
        theme_linedraw() +
        theme(
            plot.title = element_text(size = 15, hjust = 0.5, face = 'bold'),
            axis.text.x.bottom = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = 'bottom'
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

    # Return the plot
    print(p)

}
