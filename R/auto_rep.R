#' @title auto_rep
#'
#' @description Executes all GEGVIC functions to generate an HTML report with
#' the results.
#'
#' @param ge_module Logical value to determine whether this module is executed or not.
#' @param gv_module Logical value to determine whether this module is executed or not.
#' @param ic_module Logical value to determine whether this module is executed or not.
#' @param out_dir Path to determine the output directory. By default it is set
#' to the working directory.
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
#' @param kcdf 	Character string denoting the kernel to use during the non-parametric
#' estimation of the cumulative distribution function of expression levels across
#' samples when method="gsva". By default, kcdf="Poisson". It should be used when
#' input expression values are integer counts, such as those derived from RNA-seq experiments.
#' kcdf="Gaussian" which is suitable when input expression values are continuous,
#' such as microarray fluorescent units in logarithmic scale, RNA-seq log-CPMs, log-RPKMs or log-TPMs.
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
#' @param gbuild Version of the genome to work with. It can be one of the following:
#' - ‘BSgenome.Hsapiens.UCSC.hg19’
#' - ‘BSgenome.Hsapiens.UCSC.hg38’
#' - ‘BSgenome.Mmusculus.UCSC.mm10’
#' - ‘BSgenome.Mmusculus.UCSC.mm39’
#' @param mut_sigs Mutational signature matrices containing the frequencies of
#' all nucleotide changes per signature need to be indicated. GEGVIC contains the
#' matrices from COSMIC for single and double base substitutions. To choose one,
#' the user has to indicate ’COSMIC_v{XX}_{YY}BS_GRCh{ZZ}’ in the mut_sigs argument.
#' The XX is the version, that can be v2 or v3.2. YY indicates if mutations are
#' single (S) or double (D) base substitutions, while the ZZ is for the genome
#' assembly, either GRCh37 or GRCh38 for human data and mm9 or mm10 for mouse data.
#' @param tri.counts.method Normalization method. Needs to be set to either:
#' - 'default' – no further normalization.
#' - 'exome' – normalized by number of times each trinucleotide context is
#' observed in the exome.
#' - 'genome' – normalized by number of times each trinucleotide context is
#' observed in the genome.
#' - 'exome2genome' – multiplied by a ratio of that trinucleotide's occurence in
#' the genome to the trinucleotide's occurence in the exome.
#' - 'genome2exome' – multiplied by a ratio of that trinucleotide's occurence in
#' the exome to the trinucleotide's occurence in the genome.
#' - data frame containing user defined scaling factor – count data for each
#' trinucleotide context is multiplied by the corresponding value given in the
#' data frame.
#' @param indications Character vector of cancer type codes for each sample in
#' the tpm matrix.This is used by TIMER method. Indications supported
#' can be checked using immunedeconv::timer_available_cancers. Default value is
#' NULL.
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
#' @import rlang
#' @import DESeq2
#' @import SummarizedExperiment
#' @import dplyr
#' @import ggplot2
#' @import tibble
#' @import apeglm
#' @import ggrepel
#' @import clusterProfiler
#' @import GSEAmining
#' @import GSEABase
#' @import GSVA
#' @import pheatmap
#' @importFrom maftools read.maf
#' @importFrom maftools plotmafSummary
#' @importFrom maftools oncoplot
#' @import ggpubr
#' @import tidyr
#' @import deconstructSigs
#' @import ggplotify
#' @import immunedeconv
#' @import patchwork
#' @importFrom gridExtra marrangeGrob
#'
#' @examples
#' auto_rep(ge_module = TRUE,
#'          gv_module = TRUE,
#'          ic_module = TRUE,
#'          out_dir = NULL,
#'          counts = sample_counts,
#'          genes_id = 'ensembl_gene_id',
#'          metadata = sample_metadata,
#'          response = 'MSI_status',
#'          design = 'MSI_status',
#'          colors = c('orange', 'black'),
#'          ref_level = c('MSI_status', 'MSS'),
#'          shrink = 'apeglm',
#'          biomart = ensembl_biomart_GRCh38_p13,
#'          fold_change = 2,
#'          p.adj = 0.05,
#'          gmt = 'inst/extdata/c2.cp.reactome.v7.5.1.symbols.gmt',
#'          gsea_pvalue = 0.2,
#'          gsva_gmt = 'hallmark',
#'          method = 'gsva',
#'          kcdf = 'Poisson',
#'          row.names = TRUE,
#'          col.names = TRUE,
#'          muts = sample_mutations,
#'          top_genes = 10,
#'          specific_genes = NULL,
#'          compare = 'wilcox.test',
#'          p_label = 'p.format',
#'          gbuild = 'BSgenome.Hsapiens.UCSC.hg38',
#'          mut_sigs = 'COSMIC_v2_SBS_GRCh38',
#'          tri.counts.method = 'default',
#'          indications = rep('COAD', ncol(sample_counts[-1])),
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
                     kcdf = 'Poisson',
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
                                    "out_dir" = out_dir,
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
                                    "kcdf" = kcdf,
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
