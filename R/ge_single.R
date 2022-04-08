#' @title ge_single
#'
#' @description
#'
#' @param counts Data frame that contains gene expression data as raw counts.
#' @param metadata Data frame that contains supporting variables to the data.
#' @param genes_id Name of the column that contains gene identifiers. Should be
#' one of the following:'entrezgene_id', 'ensembl_gene_id' or 'hgnc_symbol'.
#' @param response Unquoted name of the variable indicating the groups to analyse.
#' @param design Variables in the design formula in the form of: 'Var1 + Var2 + ... Var_n'.
#' @param biomart Data frame containing a biomaRt query with the following
#' attributes: ensembl_gene_id, hgnc_symbol, entrezgene_id, transcript_length,
#' refseq_mrna. In the case of mus musculus data, external_gene_name must be
#' obtained and then change the column name for hgnc_symbol. Uploaded biomaRt
#' queries in GEGVIC: 'ensembl_biomartGRCh37', ensembl_biomartGRCh38_p13' and
#' 'ensembl_biomartGRCm38_p6', 'ensembl_biomartGRCm39'.
#' @param gsva_gmt Path to the gmt file that contain the gene sets of interest. By
#' default the parameter is set to 'hallmark' which provides all HALLMARK gene
#' sets from MSigDB (version 7.5.1).
#' @param method Name of the method to perform Gene set variation analysis. The
#' options are: 'gsva', 'ssgea' or 'zscore'. Default value is 'gsva'.
#' @param colors Character vector indicating the colors of the different groups
#' to compare. Default values are two: black and orange.
#' @param row.names Logical value to determine if row-names are shown in the
#' heatmap.
#' @param col.names Logical value to determine if column-names are shown in the
#' heatmap.
#'
#' @return Returns a heatmap and the expression values in a form of a matrix.
#'
#' @export
#'
#' @import rlang
#' @import DESeq2
#' @import SummarizedExperiment
#' @import dplyr
#' @import tibble
#' @import GSEABase
#' @import GSVA
#' @import pheatmap
#'
#'
#' @examples
#' gsva.res <- ge_single(counts = input_ge_module,
#'                       metadata = metadata_ge_module,
#'                       genes_id = 'entrezgene_id',
#'                       response = Response,
#'                       design = 'Response',
#'                       biomart = ensembl_biomart_GRCh38_p13,
#'                       gsva_gmt = 'hallmark',
#'                       method = 'gsva',
#'                       colors = c('black', 'orange'),
#'                       row.names = TRUE,
#'                       col.names = TRUE)
#'
ge_single <- function(counts,
                      metadata,
                      genes_id,
                      response,
                      design,
                      biomart,
                      gsva_gmt = 'hallmark',
                      method = 'gsva',
                      colors = c('black', 'orange'),
                      row.names = TRUE,
                      col.names = TRUE) {

    # Enquote response variable
    response <- rlang::enquo(response)

    # Preprocess counts data
    counts <- preprocess_ge_counts(counts = counts,
                                   genes_id = genes_id)

    # Preprocess metadata
    metadata <- preprocess_ge_meta(metadata = metadata,
                                   counts = counts) %>%
        dplyr::select(!!response)

    # Create DESeq2Dataset object
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
                                          colData = metadata,
                                          design = formula(paste('~', design,
                                                                 collapse = " ")))
    # Generate  normalized counts
    dds <- DESeq2::estimateSizeFactors(dds)

    # Transform normalized counts for data visualization
    vsd <- DESeq2::vst(dds, blind=FALSE)

    # Get expression data
    exprs.mat <- SummarizedExperiment::assay(vsd)

    # Annotate results
    ## First evaluate if genes are annotated already as hgnc_symbol
    if (genes_id == 'hgnc_symbol'){
        # If so:
        exprs.mat.annot <- as.data.frame(exprs.mat) %>%
            # Rownames to column with the name indicated in the genes_id parameter
            tibble::rownames_to_column(genes_id) %>%
            # Filter those missing gene symbols
            dplyr::filter(hgnc_symbol != '') %>%
            # Remove duplicated genes
            dplyr::distinct(hgnc_symbol, .keep_all = TRUE) %>%
            # Return hgnc_symbol column as rownames
            tibble::column_to_rownames('hgnc_symbol')

    } else {
        # If genes are identified as entrezgene_id or ensembl_gene_id:
        exprs.mat.annot <- as.data.frame(exprs.mat) %>%
            # Rownames to column with the name indicated in the genes_id parameter
            tibble::rownames_to_column(genes_id) %>%
            # Join the data frame with the GRCh38_p13 biomaRt table stored in data
            dplyr::inner_join(x = .,
                              y = biomart %>%
                                  mutate(entrezgene_id = as.character(entrezgene_id)),
                              by = genes_id) %>%
            # From all annotation columns keep only hgnc symbol column
            dplyr::select(hgnc_symbol,
                          everything(),
                          -c(entrezgene_id, ensembl_gene_id,
                             transcript_length, refseq_mrna)) %>%
            # Filter those missing gene symbols
            dplyr::filter(hgnc_symbol != '') %>%
            # Remove duplicated genes
            dplyr::distinct(hgnc_symbol, .keep_all = TRUE) %>%
            # Return hgnc_symbol column as rownames
            tibble::column_to_rownames('hgnc_symbol')
    }

    # If samples come from mouse (by the presence of one extra column)
    if('human_ortholog' %in% colnames(exprs.mat.annot) == TRUE){
        exprs.mat.annot <- exprs.mat.annot %>%
            # Rownames to column with the name indicated in the genes_id parameter
            tibble::rownames_to_column('hgnc_symbol') %>%
            # Substitute mouse gene symbol with the human homolog symbol
            dplyr::mutate(hgnc_symbol = human_ortholog) %>%
            # Filter those missing gene symbols
            dplyr::filter(hgnc_symbol != '') %>%
            # Remove duplicated genes
            dplyr::distinct(hgnc_symbol, .keep_all = TRUE) %>%
            # Remove ortholog column
            dplyr::select(-human_ortholog) %>%
            # Return hgnc_symbol column as rownames
            tibble::column_to_rownames('hgnc_symbol')
    }

    # Read genesets
    if(gsva_gmt == 'hallmark'){
        gmt <- hallmark.gmt
    } else {
        gmt <- GSEABase::getGmt(gsva_gmt)
    }

    # Calculate GSVA/ssGSEA
    gsva_temp <- GSVA::gsva(expr = as.matrix(exprs.mat.annot),
                            gset.idx.list = gmt,
                            method = method)

    # Plot Heatmap
    ## Define response level group colors in a list
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

    # Plot heatmap
    pheatmap(gsva_temp,
             scale = 'row',
             #cutree_rows = 2,
             #cutree_cols = 2,
             color = colorRampPalette(c('blue', 'grey', 'red'))(10),
             show_rownames = row.names,
             show_colnames = col.names,
             annotation_col = metadata,
             annotation_colors = pheat.anno.color,
             annotation_names_col = FALSE,
             clustering_method = 'ward.D',
             main = paste0('Samples clustering by ', toupper(method)))

    # Return results table
    return(gsva_temp)


}
