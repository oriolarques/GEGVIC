#' @title auto_rep
#'
#' @description
#'
#' @param ge_module
#' @param gv_module
#' @param ic_module
#' @param out_dir
#' @param counts Data frame that contains gene expression data as raw counts.
#' @param genes_id Name of the column that contains gene identifiers. Should be
#' one of the following:'entrez_gene_id', 'ensemblgene_id' or 'hgnc_symbol'.
#' @param metadata Data frame that contains supporting variables to the data.
#' @param response QUOTED Name of the variable indicating the groups to analyse.
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
#' @param p.adj Numeric value to define the maximum adjusted p-value to consider
#' that a gene is differentially expressed.
#' @param gmt A data frame containg the gene sets to analyse using GSEA. This
#' object should be obtained with the read.gmt function from the clusterProfiler
#' package.
#' @param gsea_pvalue Numeric value to define the adjusted pvalue cutoff during
#' GSEA. Set to 0.2 by default.
#' @param gsva_gmt Path to the gmt file that contain the gene sets of interest. By
#' default the parameter is set to 'hallmark' which provides all HALLMARK gene
#' sets from MSigDB (version 7.5.1).
#' @param method Name of the method to perform Gene set variation analysis. The
#' options are: 'gsva', 'ssgea' or 'zscore'. Default value is 'gsva'.
#' @param row.names Logical value to determine if row-names are shown in the
#' heatmap.
#' @param col.names Logical value to determine if column-names are shown in the
#' heatmap.
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
#' - Variant_Type: Type of mutation. Can be: 'SNP' (Single nucleotide polymorphism),
#' 'DNP' (Double nucleotide polymorphism), 'INS' (Insertion), 'DEL' (Deletion).
#' - Tumor_Sample_Barcode: Sample name.
#' @param top_genes Number of genes to be analysed in the mutational summary.
#' @param specific_genes Genes that will be plotted in the oncoplot.
#' @param compare A character string indicating which method to be used for
#' comparing means. Options are 't.test' and 'wilcox.test' for two groups or
#' 'anova' and 'kruskal.test' for more groups. Default value is NULL.
#' @param p_label Character string specifying label type. Allowed values include
#' 'p.signif' (shows the significance levels), 'p.format' (shows the formatted
#' p-value).
#' @param gbuild
#' @param mut_sigs
#' @param tri.counts.method
#' @param indications
#' @param cibersort Path to the CIBERSORT.R and LM22.txt files. Default value is
#' NULL.
#' @param tumor Logical value to define if samples are tumors. If so EPIC and
#' quanTIseq use a signature matrix/procedure optimized for tumor samples.
#' Default value is TRUE.
#' @param rmgenes A character vector of gene symbols. Exclude these genes from
#' the analysis. Use this to exclude e.g. noisy genes.
#' @param scale_mrna Logical. If FALSE, disable correction for mRNA content of
#' different cell types. This is supported by methods that compute an absolute
#' score (EPIC and quanTIseq). Default value is TRUE.
#' @param expected_cell_types Limit the analysis to the cell types given in this
#' list. If the cell types present in the sample are known a priori, setting
#' this can improve results for xCell
#' (see https://github.com/grst/immunedeconv/issues/1).
#' @param points Logical value to decide if points are added to the plot.
#'
#' @return
#'
#' @export
#'
#' @examples
#' auto_rep(ge_module = TRUE,
#'          gv_module = TRUE,
#'          ic_module = TRUE,
#'          out_dir = NULL,
#'          counts = input_ge_module,
#'          genes_id = 'entrezgene_id',
#'          metadata = metadata_ge_module,
#'          response = 'Response',
#'          design = 'Response',
#'          colors = c('black', 'orange'),
#'          ref_level = c('Response', 'Non_Responders'),
#'          shrink = 'apeglm',
#'          biomart = ensembl_biomart_GRCh38_p13,
#'          fold_change = 2,
#'          p.adj = 0.05,
#'          gmt = 'c7.all.v7.2.symbols.gmt',
#'          gsea_pvalue = 0.2,
#'          gsva_gmt = 'hallmark',
#'          method = 'gsva',
#'          row.names = TRUE,
#'          col.names = TRUE,
#'          muts = input_gv_module,
#'          top_genes = 10,
#'          specific_genes = NULL,
#'          compare = 'wilcox.test',
#'          p_label = 'p.format',
#'          gbuild = 'BSgenome.Hsapiens.UCSC.hg19',
#'          mut_sigs = 'COSMIC_v2_SBS_GRCh37',
#'          tri.counts.method = 'default',
#'          indications = rep('skcm', ncol(input_ge_module[-1])),
#'          cibersort = NULL,
#'          tumor = TRUE,
#'          rmgenes = NULL,
#'          scale_mrna = TRUE,
#'          expected_cell_types = NULL,
#'          points = TRUE)
#'
auto_rep <- function(ge_module = TRUE,
                     gv_module = TRUE,
                     ic_module = TRUE,
                     out_dir = NULL,
                     counts,
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
                     gsea_pvalue = 0.2,
                     gsva_gmt = 'hallmark',
                     method = 'gsva',
                     row.names = TRUE,
                     col.names = TRUE,
                     muts,
                     top_genes = 10,
                     specific_genes = NULL,
                     compare = 'wilcox.test',
                     p_label = 'p.format',
                     gbuild,
                     mut_sigs,
                     tri.counts.method = 'default',
                     indications,
                     cibersort = NULL,
                     tumor = TRUE,
                     rmgenes = NULL,
                     scale_mrna = TRUE,
                     expected_cell_types = NULL,
                     points = TRUE) {


    input <-  paste0(system.file('extdata/', package = 'GEGVIC'),
                     '/GEGVIC_auto_report.Rmd')

    if(is.null(out_dir) == TRUE){
        out_dir <- getwd()
    } else {
        out_dir <- out_dir
    }


    rmarkdown::render(input = input,
                      output_dir = out_dir,
                      params = list("ge_module" = ge_module,
                                    "gv_module" = gv_module,
                                    "ic_module" = ic_module,
                                    "counts" = counts,
                                    "genes_id" = genes_id,
                                    "metadata" = metadata,
                                    "response" = response,
                                    "design" = design,
                                    "colors" = colors,
                                    "ref_level" = ref_level,
                                    "shrink" = shrink,
                                    "biomart" = biomart,
                                    "fold_change" = fold_change,
                                    "p.adj" = p.adj,
                                    "gmt" = gmt,
                                    "gsea_pvalue" = gsea_pvalue,
                                    "gsva_gmt" = gsva_gmt,
                                    "method" = method,
                                    "row.names" = row.names,
                                    "col.names" = col.names,
                                    "muts" = muts,
                                    "top_genes" = top_genes,
                                    "specific_genes" = specific_genes,
                                    "compare" = compare,
                                    "p_label" = p_label,
                                    "gbuild" = gbuild,
                                    "mut_sigs" = mut_sigs,
                                    "tri.counts.method" = tri.counts.method,
                                    "indications" = indications,
                                    "cibersort" = cibersort,
                                    "tumor" = tumor,
                                    "rmgenes" = rmgenes,
                                    "scale_mrna" = scale_mrna,
                                    "expected_cell_types" = expected_cell_types,
                                    "points" = points))


}
