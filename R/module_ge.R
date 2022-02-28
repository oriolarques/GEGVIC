#' @title module_ge
#'
#' @description Analyses differential gene expression from RNA-seq raw counts and
#' plots PCA, volcano plots and gene set enrichment analysis (GSEA) for the
#' desired comparisons.
#'
#'
#' @param counts Data frame that contains gene expression data as raw counts.
#' @param genes_idName of the column that contains gene identifiers. Should be
#' one of the following:'entrez_gene_id', 'ensemblgene_id' or 'hgnc_symbol'.
#' @param metadata Data frame that contains supporting variables to the data.
#' @param response Unquoted name of the variable indicating the groups to analyse.
#' @param design Variables in the design formula in the form of: 'Var1 + Var2 + ... Var_n'.
#' @param colors Character vector indicating the colors of the different groups
#' to compare. Default values are two: black and orange.
#' @param ref_level Character vector where the first element is the column name
#' where the reference level is located and a second element indicating the name
#' of level to be used as a reference when calculating differential gene
#' expression.
#' @param shrink Name of the shrinkage method to apply: "apeglm", "ashr",
#' "normal" or "none". Use none to skip shrinkage. Default value is "apeglm".
#' @param biomart Data frame containing a biomaRt query with the following
#' attributes: ensembl_gene_id, hgnc_symbol, entrezgene_id, transcript_length,
#' refseq_mrna. In the case of mus musculus data, external_gene_name must be
#' obtained and then change the column name for hgnc_symbol. Uploaded biomaRt
#' queries in GEGVIC: 'ensembl_biomartGRCh37', ensembl_biomartGRCh38_p13' and
#' 'ensembl_biomartGRCm38_p6', 'ensembl_biomartGRCm39'.
#' @param fold_change An integer to define the fold change value to consider
#' that a gene is differentially expressed.
#' @param p.adjNumeric value to define the maximum adjusted p-value to consider
#' that a gene is differentially expressed.
#' @param gmt A data frame containg the gene sets to analyse using GSEA. This
#' object should be obtained with the read.gmt function from the clusterProfiler
#' package.
#' @param gsea_pvalue Numeric value to define the adjusted pvalue cutoff during
#' GSEA. Set to 0.2 by default.
#'
#' @return Returns ggplot objects containing PCA, Volcano plot and GSEA analyses.
#'
#' @export
#'
#' @import DESeq2
#' @import dplyr
#' @import ggplot2
#' @import tibble
#' @import apeglm
#' @import ggrepel
#' @import clusterProfiler
#' @import GSEAmining
#'
#' @examples
#' module_ge(counts = input_ge_module,
#'           genes_id = 'entrezgene_id',
#'           metadata = metadata_ge_module,
#'           response = Response,
#'           design = 'Response',
#'           colors = c('black', 'orange'),
#'           ref_level = c('Response', 'Non_Responders'),
#'           shrink = 'apeglm',
#'           biomart = ensembl_biomart_GRCh38_p13,
#'           fold_change = 2,
#'           p.adj = 0.05,
#'           gmt = c7.all.v7.2.symbols.gmt,
#'           gsea_pvalue = 0.2)
#'
module_ge <- function(counts,
                      genes_id,
                      metadata,
                      response,
                      design,
                      colors = c('black', 'orange'),
                      ref_level,
                      shrink = 'apeglm',
                      biomart,
                      fold_change = 2,
                      p.adj = 0.05,
                      gmt,
                      gsea_pvalue = 0.2) {

    # Get data PCA
    print('PCA')

    ## Enquote response variable
    response <- rlang::enquo(response)

    pca <- ge_pca(counts,
                  genes_id = genes_id,
                  metadata = metadata,
                  response = !!response,
                  design = design,
                  colors = colors)

    print(pca)

    # Run differential gene expression analysis
    print('Differential gene expression analysis')

    results.dds <- ge_diff_exp(counts = counts,
                               genes_id = genes_id,
                               metadata = metadata,
                               design = design,
                               ref_level = ref_level,
                               shrink = shrink)


    # Annotate gene symbols
    print('Annotate gene symbols')

    annot.res <- ge_annot(results_dds = results.dds,
                          genes_id = genes_id,
                          biomart = biomart)

    # Create volcano plot
    print('Volcano plots')

    ge_volcano(annot_res = annot.res,
               fold_change = fold_change,
               p.adj = p.adj)

    # Obtain GSEA results
    print('GSEA')

    ge_gsea(annot_res = annot.res,
            gmt = gmt,
            gsea_pvalue = gsea_pvalue)


}

