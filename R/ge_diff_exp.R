#' @title ge_diff_exp
#'
#' @description
#'
#' @param counts Data frame that contains gene expression data as raw counts.
#' @param genes_id Name of the column that contains gene identifiers. Should be
#' one of the following:'entrez_gene_id', 'ensemblgene_id' or 'hgnc_symbol'.
#' @param metadata Data frame that contains supporting variables to the data.
#' @param design Variables in the design formula in the form of: 'Var1 + Var2 + ... Var_n'.
#' @param ref_level Character vector where the first element is the column name
#' where the reference level is located and a second element indicating the name
#' of level to be used as a reference when calculating differential gene
#' expression.
#'
#' @return Returns a list of differential gene expression results, one for each
#' level comparison, in form of DESeqResults objects.
#'
#' @export
#'
#' @import DESeq2
#' @import apeglm
#' @import dplyr
#' @import ggplot2
#' @import tibble
#'
#' @examples
#' results_dds <- ge_diff_exp(counts = input_ge_module,
#'                            genes_id = 'entrezgene_id',
#'                            metadata = metadata_ge_module,
#'                            design = 'Response',
#'                            ref_level = c('Response', 'Non_Responders'))
#'
ge_diff_exp <- function(counts,
                        genes_id,
                        metadata,
                        design,
                        ref_level) {

    # Input preprocessing -----------------------------------------------------
    # Preprocess counts data
    counts <- preprocess_ge_counts(counts = counts,
                                   genes_id = genes_id)

    # Preprocess meta data
    metadata <- preprocess_ge_meta(metadata = metadata)

    # Create DESeq2Dataset object ---------------------------------------------
    dds <- DESeqDataSetFromMatrix(countData = counts,
                                  colData = metadata,
                                  design = formula(paste('~', design, collapse = " ")))

    # Define the reference level in the grouping variable
    dds[[ref_level[1]]] <- relevel(dds[[ref_level[1]]], ref = ref_level[2])

    # Compute differential gene expression using DESeq ------------------------
    dds <- DESeq(dds)

    # Shrinkage estimation ----------------------------------------------------
    ## Create an empty list to store the different comparisons in the case
    ## there are multiple levels in the variable of study
    results_dds <- list()

    # Iterate over the different comparisons
    for (i in seq(from = 2, to = length(resultsNames(dds)))) {
        # Calculate shrinkage estimators and save results in the results_dds list
        results_dds[[i-1]] <- lfcShrink(dds = dds,
                                        coef = resultsNames(dds)[i],
                                        type =  'apeglm')
        # Add the name of each comparison to the corresponding element
        names(results_dds)[i-1] <- resultsNames(dds)[i]
    }

    return(results_dds)
}
