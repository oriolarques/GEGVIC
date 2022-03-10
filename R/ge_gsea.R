#' @title ge_gsea
#'
#' @description
#'
#' @param annot_res The output of ge_annot. List containing data frames of
#' differential gene expression results between the different groups.
#' @param gmt Path to the gmt file that contain the gene sets of interest.
#' @param gsea_pvalue Numeric value to define the adjusted pvalue cutoff during
#' GSEA. Set to 0.2 by default.
#'
#' @return Returns ggplot objects and a list of gseaResult objects.
#'
#' @export
#'
#' @import dplyr
#' @import clusterProfiler
#' @import GSEAmining
#'
#' @examples
#' results_dds <- ge_diff_exp(counts = input_ge_module,
#'                            genes_id = 'entrezgene_id',
#'                            metadata = metadata_ge_module,
#'                            design = 'Response',
#'                            ref_level = c('Response', 'Non_Responders'),
#'                            shrink = 'none')
#' annot.res <- ge_annot(results_dds = results_dds,
#'                       genes_id = 'entrezgene_id',
#'                       biomart = ensembl_biomart_GRCh38_p13)
#' gsea.res <- ge_gsea(annot_res = annot.res,
#'                     gmt = 'inst/extdata/c7.all.v7.2.symbols.gmt',
#'                     gsea_pvalue = 0.2)
#'
ge_gsea <- function(annot_res,
                    gmt,
                    gsea_pvalue = 0.2) {

    # 0. If samples come from mouse (by the presence of one extra column)
    if('human_ortholog' %in% colnames(annot_res[[1]]) == TRUE){
        # Iterate over annotated results list
        for (i in seq_along(annot_res)){
            # By each data frame in annot_res
            annot_res[[i]] <- annot_res[[i]] %>%
                # Substitute mouse gene symbol with the human homolog symbol
                dplyr::mutate(hgnc_symbol = human_ortholog) %>%
                # Filter those missing gene symbols
                dplyr::filter(hgnc_symbol != '') %>%
                # Remove duplicated genes
                dplyr::distinct(hgnc_symbol, .keep_all = TRUE)
        }

    }

    # 1. Obtain geneLists
        ## Create a list to store as many geneList as conditions
        geneLists <- list()

        ## Create two temporal objects
        temp_df <- NULL
        temp_gs <- NULL

        ## Create geneLists
        ### Iterate over the annotated results list
        for (i in seq_along(annot_res)) {
            # Get the data frame with gene expression and sort by fold change
            temp_df <- annot_res[[i]] %>%
                #Filter genes whose log2FoldChange == NA
                filter(!is.na(.data$log2FoldChange)) %>%
                arrange(desc(.data$log2FoldChange))

            # Create a vector with the log2FoldChange
            temp_gs <- temp_df$log2FoldChange
            # Get the gene symbols of each log2FoldChange
            names(temp_gs) <- as.character(temp_df$hgnc_symbol)

            # Store the geneList
            geneLists[[i]] <- temp_gs
            names(geneLists)[i] <- paste0('geneList_', names(annot_res)[i])
        }

        # 2. Perform GSEA
        ## Create an empty list to store GSEA results
        GSEA.res <- list()
        temp_gsea <- NULL

        ## Read the gmt file
        gmt <- clusterProfiler::read.gmt(gmt)

        ## Execute GSEA
        ### Iterate over the annotated results list
        for (i in seq_along(annot_res)) {
            temp_gsea <- clusterProfiler::GSEA(geneList = geneLists[[i]],
                                               TERM2GENE = gmt,
                                               pvalueCutoff = gsea_pvalue)
            # Delay the next process 0.5 seconds
            Sys.sleep(0.5)
            # Save the GSEA results
            GSEA.res[[i]] <- temp_gsea
            names(GSEA.res)[i] <- paste0('GSEA_', names(annot_res)[i])

            # Check if GSEA result is empty and print a message
            if (nrow(GSEA.res[[i]]@result) == 0) {
                print('No gene sets are enriched under specific pvalueCutoff')
            } else {
                # 2.1. GSEAmining
                # Filter gene sets to analyse the top ones
                gs.filt <- GSEA.res[[i]]@result %>%
                    dplyr::arrange(desc(.data$NES)) %>%
                    dplyr::mutate(group = ifelse(test = .data$NES > 0,
                                                 yes = 'Positive',
                                                 no = 'Negative')) %>%
                    dplyr::group_by(.data$group) %>%
                    dplyr::filter(.data$p.adjust < gsea_pvalue) %>%
                    dplyr::top_n(., n = 20, wt = abs(.data$NES)) %>%
                    dplyr::ungroup(.) %>%
                    dplyr::select(.data$ID,
                                  .data$NES,
                                  .data$p.adjust,
                                  .data$leading_edge,
                                  .data$core_enrichment)

                # 2.2. Bubble plot
                print(
                    gs.filt %>%
                        separate(leading_edge , into= 'tags', sep=',') %>%
                        separate(tags, into = c('tags', 'core_perc'), sep='=') %>%
                        separate(core_perc, into = 'core_perc', sep='%') %>%
                        mutate(core_perc = as.numeric(core_perc)) %>%
                        dplyr::select(-tags) %>%
                        ggplot(.,
                               aes(x= NES,
                                   y=reorder(ID, NES),
                                   size= core_perc,
                                   colour = p.adjust))+
                        geom_point(alpha=0.5)+
                        geom_vline(xintercept = 0)+
                        scale_color_gradient(low = "#FF9900", high = "#FF3300")+
                        labs(title =  names(annot_res)[i],
                             size='% of genes in\n leading edge', colour = 'p.adjust')+
                        theme_bw()+
                        theme(panel.grid = element_blank(),
                              axis.text = element_text(size=12, face = "bold"),
                              axis.title.y = element_blank(),
                              axis.title.x = element_text(size=15),
                              axis.text.x = element_text(size=10),
                              legend.title = element_text(face='bold', size =8),
                              legend.text = element_text(size =7))
                )

                # 2.3. Cluster gene sets
                gs.cl <- GSEAmining::gm_clust(df = gs.filt)

                # Plot cluster
                GSEAmining::gm_dendplot(df = gs.filt,
                                        hc = gs.cl)

                # 2.4. Plot enriched terms in gene sets names
                print(GSEAmining::gm_enrichterms(df = gs.filt,
                                                 hc = gs.cl))

                # 2.5. Plot enriched cores (leading edge analysis)
                print(GSEAmining::gm_enrichcores(df = gs.filt,
                                                 hc = gs.cl))

                }
        }

        return(GSEA.res)

}
