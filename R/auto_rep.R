#' auto_rep
#'
#' @param counts
#' @param genes_id
#' @param metadata
#' @param muts
#' @param response
#' @param design
#' @param colors
#' @param ref_level
#' @param shrink
#' @param biomart
#' @param fold_change
#' @param p.adj
#' @param gmt
#' @param gsea_pvalue
#' @param gsva_gmt
#' @param method
#' @param row.names
#' @param col.names
#' @param top_genes
#' @param specific_genes
#' @param compare
#' @param p_label
#' @param gbuild
#' @param mut_sigs
#' @param tri.counts.method
#' @param indications
#' @param cibersort
#' @param tumor
#' @param rmgenes
#' @param scale_mrna
#' @param expected_cell_types
#'
#' @return
#' @export
#'
#' @examples
#' prova(counts = input_ge_module,
#'       genes_id = 'entrezgene_id',
#'       metadata = met,
#'       response = 'Response',
#'       design = 'Response',
#'       colors = c('black', 'orange', 'purple'),
#'       ref_level = c('Response', 'Non_Responders'),
#'       shrink = 'none',
#'       biomart = ensembl_biomart_GRCh38_p13,
#'       fold_change = 2,
#'       p.adj = 0.05,
#'       gmt = 'c7.all.v7.2.symbols.gmt',
#'       gsea_pvalue = 0.2,
#'       gsva_gmt = 'hallmark',
#'       method = 'gsva',
#'       row.names = TRUE,
#'       col.names = TRUE,
#'       muts = input_gv_module,
#'       top_genes = 10,
#'       specific_genes = NULL,
#'       compare = 'wilcox.test',
#'       p_label = 'p.format',
#'       gbuild = 'BSgenome.Hsapiens.UCSC.hg19',
#'       mut_sigs = 'COSMIC_v2_SBS_GRCh37',
#'       tri.counts.method = 'default',
#'       indications = rep('skcm', ncol(input_ge_module[-1])),
#'       cibersort = NULL,
#'       tumor = TRUE,
#'       rmgenes = NULL,
#'       scale_mrna = TRUE,
#'       expected_cell_types = NULL)
#'
auto_rep <- function(counts,
                     genes_id,
                     metadata,
                     muts,
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
                     expected_cell_types = NULL) {


    input <-  paste0(system.file('extdata/', package = 'GEGVIC'),
                     '/GEGVIC_auto_report.Rmd')

    out_dir <- getwd()

    rmarkdown::render(input = input,
                      output_dir = out_dir,
                      params = list("counts" = counts,
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
                                    "expected_cell_types" = expected_cell_types))


}
