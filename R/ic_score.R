#' @title ic_score
#'
#' @description Calculates and plots immunophenogram (IPG) and immunophenoscores
#' (IPS) for each sample and each group of study.
#'
#' @param tpm Output from the ic_raw_to_tpm function. This is a matrix
#' containing expression counts as TPM with HGNC gene symbols as rownames
#' and samples identifiers as colnames.
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
#' @return Returns ggplot objects and a data frame.
#'
#' @export
#'
#' @import dplyr
#' @import ggplot2
#' @import rlang
#' @import grid
#' @import ggpubr
#' @import patchwork
#' @import ggplotify
#' @importFrom gridExtra marrangeGrob
#'
#' @examples
#' tpm <- ic_raw_to_tpm(counts = input_ge_module,
#'                      genes_id = 'entrezgene_id',
#'                      biomart = ensembl_biomart_GRCh38_p13)
#' ips <- ic_score(tpm = tpm,
#'                 metadata = sample_metadata,
#'                 response = MSI_status,
#'                 compare = 'wilcox.test',
#'                 p_label = 'p.format',
#'                 colors = c('orange', 'black'))
#'
ic_score <- function(tpm,
                     metadata,
                     response,
                     compare = NULL,
                     p_label = 'p.format',
                     colors = c('orange', 'black')) {

    # Create a list to store phenograms
    ipheno_list <- list()
    # Create a list to store the results
    ic.score.results <- list()

    # Preliminary functions ---------------------------------------------------
    ## To calculate Immunophenoscore
    ipsmap<- function (x) {
        if (x<=0) {
            ips<-0
        } else {
            if (x>=3) {
                ips<-10
            } else {
                ips<-round(x*10/3, digits=0)
            }
        }
        return(ips)
    }

    ## To assign colors
    my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 1000)
    mapcolors<-function (x) {
        za<-NULL
        if (x>=3) {
            za=1000
        } else {
            if (x<=-3) {
                za=1
            } else {
                za=round(166.5*x+500.5,digits=0)
            }
        }
        return(my_palette[za])
    }
    my_palette2 <- colorRampPalette(c("black", "white"))(n = 1000)
    mapbw<-function (x) {
        za2<-NULL
        if (x>=2) {
            za2=1000
        } else {
            if (x<=-2) {
                za2=1
            } else {
                za2=round(249.75*x+500.5,digits=0)
            }
        }
        return(my_palette2[za2])
    }

    # Load and preprocess necessary data --------------------------------------

    ## Transform TPM data into log2(TPM+1)
    gene_expression <- as.data.frame(log2(tpm + 1))
    sample_names <- names(gene_expression)

    ## Read IPS genes and corresponding weights from tab-delimited text file "IPS_genes.txt"
    IPSG <- GEGVIC:::IPS_genes
    unique_ips_genes <- as.vector(unique(IPSG$NAME))

    ## Create variables
    IPS <- NULL
    MHC <- NULL
    CP <- NULL
    EC <- NULL
    SC <- NULL
    AZ <- NULL

    # Detect missing genes ----------------------------------------------------
    ## Gene names in expression file
    #GVEC <- row.names(gene_expression)
    ## Genes names in IPS genes file
    #VEC <- as.vector(IPSG$GENE)
    ## Match IPS genes with genes in expression file
    #ind <- which(is.na(match(VEC,GVEC)))
    ## List genes missing or differently named
    #MISSING_GENES <- VEC[ind]
    #dat <- IPSG[ind,]
    #if (length(MISSING_GENES)>0) {
    #    cat("differently named or missing genes: ",MISSING_GENES,"\n")
    #}
    #for (x in 1:length(ind)) {
    #    print(IPSG[ind,])
    #}

    ## Enquote response variable
    response <- enquo(response)


    # Calculate ImmunoPhenoGram -----------------------------------------------
    for (i in 1:length(sample_names)) {
        GE <- gene_expression[[i]]
        mGE <- mean(GE)
        sGE <- sd(GE)
        Z1 <- (gene_expression[as.vector(IPSG$GENE),i]-mGE)/sGE
        W1 <- IPSG$WEIGHT
        WEIGHT <- NULL
        MIG <- NULL
        k <- 1
        for (gen in unique_ips_genes) {
            MIG[k] <- mean(Z1[which (as.vector(IPSG$NAME)==gen)],na.rm=TRUE)
            WEIGHT[k] <- mean(W1[which (as.vector(IPSG$NAME)==gen)])
            k <- k+1
        }

        ## Chunk added to avoid problems with data coming from mouse
        ## such as genes missing in the process of finding orthologs
        # In case there is an element as NA due to missing genes in data
        for(x in seq_along(MIG)){
            if(is.na(MIG[x])){
                MIG[x]<-0
            }
        }


        WG<-MIG*WEIGHT
        MHC[i] <- mean(WG[1:10])
        CP[i] <- mean(WG[11:20])
        EC[i] <- mean(WG[21:24])
        SC[i] <- mean(WG[25:26])
        AZ[i] <- sum(MHC[i],CP[i],EC[i],SC[i])
        IPS[i] <- ipsmap(AZ[i])

        # Plot Immunophenogram  -----------------------------------------------
        data_a <-  data.frame(start = c(0,2.5,5,7.5,10,15,seq(20,39),0,10,20,30),
                              end = c(2.5,5,7.5,10,15,seq(20,40),10,20,30,40),
                              y1=c(rep(2.6,26),rep(0.4,4)),
                              y2=c(rep(5.6,26),rep(2.2,4)),
                              z=c(MIG[c(21:26,11:20,1:10)],
                                  EC[i],SC[i],CP[i],MHC[i]),
                              vcol=c(unlist(lapply(MIG[c(21:26,11:20,1:10)],
                                                   mapcolors)),
                                     unlist(lapply(c(EC[i],SC[i],CP[i],MHC[i]),
                                                   mapbw))),
                              label = c(unique_ips_genes[c(21:26,11:20,1:10)],
                                        "EC","SC","CP","MHC"))

        data_a$label <- factor(data_a$label,
                               levels=unique(data_a$label))

        ### Plot IPG per sample
        plot_a <- ggplot() +
            geom_rect(data=data_a,
                      mapping=aes(xmin=start, xmax=end, ymin=y1, ymax=y2, fill=label),
                      size=0.5,color="black", alpha=1) +
            coord_polar() +
            scale_y_continuous(limits = c(0, 6)) +
            scale_fill_manual(values = as.vector(data_a$vcol),guide='none') +
            theme_bw() +
            theme(panel.spacing = unit(0, 'mm'),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(colour = "white"),
                  axis.text=element_blank(),
                  axis.ticks= element_blank())

        # Get the Response name for the sample
        samp.resp <- metadata %>%
            filter(Samples == sample_names[i]) %>%
            pull(!!response)

        ### Add Sample name as plot title
        plot_a <- plot_a  +
            labs(title = sample_names[i],
                 subtitle = samp.resp) +
            theme(plot.title = element_text(size = 15, vjust = -1, hjust = 0.5),
                  plot.subtitle = element_text(size = 15, vjust = -1, hjust = 0.5),
                  plot.margin=unit(c(0,0,0,0),"mm"))

        ### Plot Legend sample-wise (averaged) z-scores -----------------------
        data_b <- data.frame (start = rep(0,23),
                              end = rep(0.7,23),
                              y1=seq(0,22,by=1),
                              y2=seq(1,23,by=1),
                              z=seq(-3,3,by=6/22),
                              vcol=c(unlist(lapply(seq(-3,3,by=6/22),mapcolors))),
                              label = LETTERS[1:23])
        data_b_ticks <- data.frame(x = rep(1.2, 7),
                                   value = seq(-3,3, by=1),
                                   y = seq(0,6, by=1)*(22/6) +0.5)
        legendtheme <- theme(plot.margin = unit(c(2,0,2,0),"inch"),
                             panel.spacing = unit(0,"null"),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.border = element_blank(),
                             panel.background = element_blank(),
                             axis.line = element_line(colour = "white"),
                             axis.text=element_blank(),
                             axis.ticks= element_blank(),
                             axis.title.x=element_blank())

        plot_b <- ggplot(hjust=0) +
            geom_rect(data=data_b,
                      mapping=aes(xmin=start, xmax=end, ymin=y1, ymax=y2, fill=label),
                      size=0.5,color="black", alpha=1) +
            scale_x_continuous(limits = c(0, 1.5),
                               expand = c(0,0)) +
            scale_fill_manual(values =as.vector(data_b$vcol),
                              guide='none') +
            geom_text(data=data_b_ticks,
                      aes(x=x, y=y, label=value),
                      hjust="inward",
                      size=4) +
            theme_bw() +
            legendtheme +
            ylab("Sample-wise (averaged) z-score")

        ### Plot Legend weighted z-scores -------------------------------------
        data_c <- data.frame (start = rep(0,23),
                              end = rep(0.7,23),
                              y1=seq(0,22,by=1),
                              y2=seq(1,23,by=1),
                              z=seq(-2,2,by=4/22),
                              vcol=c(unlist(lapply(seq(-2,2,by=4/22),mapbw))),
                              label = LETTERS[1:23])
        data_c_ticks <- data.frame(x = rep(1.2, 5),
                                   value = seq(-2,2, by=1),
                                   y = seq(0,4, by=1)*(22/4) +0.5)
        plot_c<-ggplot() +
            geom_rect(data=data_c,
                      mapping=aes(xmin=start, xmax=end, ymin=y1, ymax=y2, fill=label),
                      size=0.5,color="black",
                      alpha=1) +
            scale_x_continuous(limits = c(0, 1.5),
                               expand = c(0,0)) +
            scale_fill_manual(values =as.vector(data_c$vcol),
                              guide='none') +
            geom_text(data=data_c_ticks,
                      aes(x=x, y=y, label=value),
                      hjust="inward",
                      size=4) +
            theme_bw() +
            legendtheme +
            ylab("Weighted z-score")

        ## Save plot to file (1 pdf file for each sample)
        #file_name<-paste("IPS_",sample_names[i],".pdf",sep="")
        #pdf(file_name, width=10, height=8)
        #grid.arrange(plot_a,plot_b,plot_c, ncol=3, widths=c(0.8,0.1,0.1))
        #dev.off()

        ## Save each patchwork IPG plot for each patient in a list
        ipheno_list[[i]] <- (plot_a + plot_b + plot_c) +
            patchwork::plot_layout(widths = c(1,0.5,0.5))
        ## Indicate the name of the patient in the list
        names(ipheno_list)[i] <- sample_names[i]

    }

    # Save a pdf report with each IPG per sample
    report <- lapply(ipheno_list, function(x) ggplotify::as.ggplot(plot = x,
                                                                   scale = 1.1))
    report <- marrangeGrob(grobs = report, ncol = 1, nrow = 1)

    ggsave('immunophenogram_report.pdf', report)

    # # Plot immunophenogram by group -------------------------------------------
    #
    # ## Enquote response variable
    # response <- enquo(response)
    # ## Get the levels of response (grouping) variable
    # levels.resp <- levels(metadata %>%
    #                           dplyr::mutate_all(as.factor) %>%
    #                           dplyr::select(!!response) %>%
    #                           pull(.))
    # ## Create a list to store the plots for each level in response variable
    # list.resp <- list()
    #
    # ## Iterate over levels in response variable
    # for(i in seq_along(levels.resp)){
    #     # Get the sample names that belong to one of the levels in response variable
    #     temp_samp.resp <- metadata %>%
    #         dplyr::filter(!!response == levels.resp[i]) %>%
    #         dplyr::select(Samples) %>% pull(.)
    #     # Arrange IPG for all samples in one level of response variable
    #     list.resp[[i]] <- ggpubr::ggarrange(plotlist = ipheno_list[temp_samp.resp],
    #                                         labels = levels.resp[i],
    #                                         hjust = -1,
    #                                         vjust = 0.5)
    #
    #
    # }

    # Plot immunophenoscore by group ------------------------------------------
    ## Create a data frame to store
    DF<-data.frame(Samples=sample_names, MHC=MHC, EC=EC, SC=SC,
                   CP=CP, AZ=AZ, IPS=IPS)

    ## Add metadata information
    DF <- DF %>%
        dplyr::left_join(x = .,
                         y = metadata,
                         by = 'Samples') %>%
        dplyr::group_by(!!response)


    #write.table(DF,file="IPS.txt",row.names=FALSE, quote=FALSE,sep="\t")

    ## Plot IPS and subplots
    ### Define the titles of the subplots
    plot_names <- c('MHC' = 'MHC: MHC molecules',
                    'CP' = 'CP: Immunomodulators',
                    'EC' = 'EC: Effector cells',
                    'SC' = 'SC: Suppressor cells')
    ### Subplots
    subplots <- DF %>%
        # Reshape the data to long format
        dplyr::select(Samples, MHC, EC, SC, CP, !!response) %>% # Do not include AZ or IPS
        tidyr::pivot_longer(cols = -c(Samples, !!response),
                            names_to = 'Cell_type',
                            values_to = 'Score') %>%
        # Re-order the levels
        dplyr::mutate(Cell_type = factor(Cell_type, levels = c('EC', 'SC',
                                                               'CP', 'MHC'))) %>%
        dplyr::group_by(Cell_type) %>%
        # Plot
        ggplot(., aes(x = !!response,
                      y = Score,
                      col = !!response)) +
        # Geometric objects
        geom_violin() +
        geom_boxplot(width = 0.1, outlier.shape = NA) +
        geom_point(alpha = 0.5, position = position_jitter(0.2)) +

        # Define colors
        scale_color_manual(values = colors) +

        # Themes
        theme_bw() +
        theme(text = element_text(size = 15),
              plot.title = element_text(size = 15, hjust = 0.5, face = 'bold'),
              axis.text.x.bottom = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = 'bottom',
              strip.background = element_rect(
                  color="black", fill="black", size=1.5, linetype="solid"),
              strip.text = element_text(color = 'white')) +
        facet_wrap(~ Cell_type,
                   scales = 'free_y',
                   ncol = 2,
                   labeller = as_labeller(plot_names))

    ### IPS plot
    ips.plot <- ggplot(DF, aes(x = !!response,
                               y = IPS,
                               col = !!response)) +
        # Geometric objects
        geom_violin() +
        geom_boxplot(width = 0.1, outlier.shape = NA) +
        geom_point(alpha = 0.5, position = position_jitter(0.2)) +

        # Define colors
        scale_color_manual(values = colors) +

        # Title
        ggtitle('ImmunophenoScore') +

        # Themes
        theme_bw() +
        theme(text = element_text(size = 15),
              plot.title = element_text(size = 15, hjust = 0.5, face = 'bold'),
              axis.text.x.bottom = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = 'bottom'
        )
    if(is.null(compare) == FALSE){
        subplots <- subplots +
            ggpubr::stat_compare_means(method = compare,
                                       label = p_label,
                                       #size = 5,
                                       label.y.npc = 1,
                                       label.x.npc = 0.5,
                                       show.legend = FALSE)
        ips.plot <- ips.plot +
            ggpubr::stat_compare_means(method = compare,
                                       label = p_label,
                                       size = 5,
                                       label.y.npc = 1,
                                       label.x.npc = 0.5,
                                       show.legend = FALSE)

    }

    ### Add subplots and IPS plot together
    final.ips.plot <- ggpubr::ggarrange(plotlist = list(subplots, ips.plot),
                                        hjust = -1,
                                        vjust = 0.5,
                                        common.legend	= TRUE,
                                        legend = 'bottom')


    # Plot IPS results ------------------------------------------------

    ## Add together IPG from all levels in response variable
    # ipg.plots <- patchwork::wrap_plots(list.resp)
    # ipg.legends <- patchwork::wrap_plots(plot_b, plot_c)
    # layout <- c(
    #     patchwork::area(t = 1, l = 1, b = 5, r = 4),
    #     patchwork::area(t = 2, l = 5, b = 4, r = 5)
    # )
    #
    # ic.score.results[[1]] <- patchwork::wrap_plots(ipg.plots, ipg.legends,
    #                                                design = layout)
    #

    ## Get IPS comparison between samples in different levels of response variable
    # ic.score.results[[2]] <- final.ips.plot
    # names(ic.score.results) <- c('immunophenoGram', 'immunophenoScore')
    #print(ic.score.results)

    print(final.ips.plot)

    # Return results data frame
    return(DF)
}
