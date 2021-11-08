#' @title ge_volcano
#'
#' @description Plots volcano plots from the outputs of the ge_annot() function.
#'
#' @param annot_res The output of ge_annot. List containing data frames of
#' differential gene expression results between the different groups.
#' @param fold_change An integer to define the fold change value to consider
#' that a gene is differentially expressed.
#' @param p.adj Numeric value to define the maximum adjusted p-value to consider
#' that a gene is differentially expressed.
#'
#' @return Returns ggplot graphs.
#'
#' @export
#'
#' @import dplyr
#' @import ggplot2
#' @import ggrepel
#'
#' @examples
#' results_dds <- ge_diff_exp(counts = input_ge_module,
#'                            genes_id = 'entrezgene_id',
#'                            metadata = metadata_ge_module,
#'                            design = 'Response',
#'                            ref_level = c('Response', 'Non_Responders'),
#'                            shrink = 'apeglm')
#' annot.res <- ge_annot(results_dds = results_dds, genes_id = 'entrezgene_id')
#' ge_volcano(annot_res = annot.res, fold_change = 2, p.adj = 0.05)

ge_volcano <- function(annot_res,
                       fold_change = 2,
                       p.adj = 0.05) {

    # Create a temporal object to store the plots within the loop
    temp_plot <- NULL

    # Iterate over the annotated results list
    for (i in seq_along(annot_res)) {

        # Arrange the data
        volcano.plot <- annot_res[[i]] %>%
            # filter padj == NA
            filter(!is.na(.data$padj)) %>%
            # Arrange data in descending order
            arrange(desc(.data$log2FoldChange))

        # Make a volcano plot
        temp_plot <- ggplot(volcano.plot, aes(x = .data$log2FoldChange,
                                       y = -log10(.data$padj))) +
            # Plot all the genes in grey
            geom_point(alpha = 0.55, col = 'grey')+
            # Plot the significantly up-regulated genes
            geom_point(data = subset(x = volcano.plot,
                                     subset = log2FoldChange > log2(fold_change) & padj < p.adj),
                       mapping = aes(x = .data$log2FoldChange,
                                     y = -log10(.data$padj)),
                       col = 'red3',
                       alpha = 0.6) +
            # Plot the significantly down-regulated genes
            geom_point(data = subset(x = volcano.plot,
                                     subset = log2FoldChange < -log2(fold_change) & padj < p.adj),
                       mapping = aes(x = .data$log2FoldChange,
                                     y = -log10(.data$padj)),
                       col = 'blue3',
                       alpha = 0.6) +
            # Label the top 10 up-regulated genes
            geom_text_repel(data = top_n(subset(x = volcano.plot,
                                                subset = log2FoldChange > log2(fold_change) &
                                                    padj < p.adj),
                                         n = - 10,
                                         wt = .data$padj),
                            mapping = aes(label = .data$hgnc_symbol),
                            col = 'darkred',
                            box.padding = 0.5,
                            max.overlaps = 1000) +
            # Label the top 10 down-regulated genes
            geom_text_repel(data = top_n(subset(x = volcano.plot,
                                                subset = log2FoldChange < -log2(fold_change) &
                                                    padj < p.adj),
                                         n = - 10,
                                         wt = .data$padj),
                            mapping = aes(label = .data$hgnc_symbol),
                            col = 'darkblue',
                            box.padding = 0.5,
                            max.overlaps = 1000) +

            # Plot lines to identify the log2FoldChange -2 and 2 and the -log10(padj) < 0.01
            geom_vline(xintercept = c(-log2(fold_change),
                                      log2(fold_change)),
                       linetype = 'dashed')+
            geom_hline(yintercept = -log10(p.adj),
                       linetype = 'dashed')+

            # Define dynamic breaks for every log2 unit in fold change
            scale_x_continuous(breaks = seq(from = round(min(volcano.plot$log2FoldChange),0),
                                            to = round(max(volcano.plot$log2FoldChange), 0),
                                            by = 1)) +

            # Define the title, theme and axes
            ggtitle(names(annot_res)[i])+
            xlab('Log2 Fold Change')+
            ylab('-log10(Adjusted p-value)')+
            theme_bw()+
            theme(#legend.position = 'none',
                  panel.border = element_rect(fill=NA, size=2),
                  panel.grid = element_blank(),
                  axis.text = element_text(size=15, face = "bold"),
                  axis.title = element_text(size=15, face = 'bold'),
                  plot.title = element_text(size=15, face = 'bold', hjust = 0.5))

        # Plot the results
        print(temp_plot)

    }

}
