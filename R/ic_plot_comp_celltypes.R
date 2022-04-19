#' @title ic_plot_comp_celltypes
#'
#' @description Plots immune cell fractions to be compared within each sample.
#' Results from three different methods (CIBERSORT, EPIC and QUANTISEQ) are
#' presented by each group of samples.
#'
#' @param df Output ot ic_deconv function. A data frame containing the predictions
#' of the six methods included in the immunedeconv package. First column must
#' indicate and be named as cell_type, method and the rest being the numeric
#' prediction for each sample.
#' @param metadata Data frame that contains supporting variables to the data.
#' @param response Unquoted name of the variable indicating the groups to analyse.
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
#' ic_plot_comp_celltypes(df = ic.pred,
#'                      metadata = sample_metadata,
#'                      response = MSI_status)
#'
ic_plot_comp_celltypes <- function(df,
                                   metadata,
                                   response) {

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
        # Select methods that allow within samples comparison
        dplyr::filter(method %in% c('EPIC', 'QUANTISEQ', 'CIBERSORT')) %>%
        # Group by Sample, cell_type and method
        dplyr::group_by(Samples, cell_type, method) %>%
        # Sum the prediction values within different cell_types of the same
        ## category in each patient
        dplyr::summarise(cell_type = cell_type,
                         response = (!!response),
                         estimation = sum(estimation)) %>%
        # Eliminate duplicated rows
        dplyr::distinct()

    # Plot Within sample comparison of cell types
    p <- ggplot(temp_df, aes(x = Samples,
                             y = estimation,
                             fill = cell_type)) +

        # Geometric objects
        geom_bar(stat = 'identity') +

        # Define fill colors using the Accent palette from ggpubr package
        scale_fill_manual(values = ggpubr::get_palette(palette = 'Accent',
                                                       k = length(unique(
                                                           temp_df$cell_type
                                                           )))) +
        # Expand columns to fill margins
        scale_y_continuous(expand = c(0,0)) +

        # Title
        ggtitle('Immune composition:\nWithin sample comparison of cell types') +

        # Themes
        theme_bw() +
        theme(
            plot.title = element_text(size = 15, hjust = 0.5, face = 'bold'),
            axis.text.x.bottom = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = 5, angle = 45, hjust = 1),
            strip.background = element_rect(
                color="black", fill="black", size=1.5, linetype="solid"),
            strip.text = element_text(color = 'white')
        ) +

        # Faceting
        facet_grid(method ~ response,
                   scales = 'free_x',
                   switch = 'y')

    # Return the plot
    print(p)



}
