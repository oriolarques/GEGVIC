% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ge_single.R
\name{ge_single}
\alias{ge_single}
\title{ge_single}
\usage{
ge_single(
  counts,
  metadata,
  genes_id,
  response,
  design,
  biomart,
  gsva_gmt = "hallmark",
  method = "gsva",
  kcdf = "Gaussian",
  colors = c("orange", "black"),
  row.names = TRUE,
  col.names = TRUE
)
}
\arguments{
\item{counts}{Data frame that contains gene expression data as raw counts.}

\item{metadata}{Data frame that contains supporting variables to the data.}

\item{genes_id}{Name of the column that contains gene identifiers. Should be
one of the following:'entrezgene_id', 'ensembl_gene_id' or 'hgnc_symbol'.}

\item{response}{Unquoted name of the variable indicating the groups to analyse.}

\item{design}{Variables in the design formula in the form of: 'Var1 + Var2 + ... Var_n'.}

\item{biomart}{Data frame containing a biomaRt query with the following
attributes: ensembl_gene_id, hgnc_symbol, entrezgene_id, transcript_length,
refseq_mrna. In the case of mus musculus data, external_gene_name must be
obtained and then change the column name for hgnc_symbol. Uploaded biomaRt
queries in GEGVIC: 'ensembl_biomartGRCh37', ensembl_biomartGRCh38_p13' and
'ensembl_biomartGRCm38_p6', 'ensembl_biomartGRCm39'.}

\item{gsva_gmt}{Path to the gmt file that contain the gene sets of interest. By
default the parameter is set to 'hallmark' which provides all HALLMARK gene
sets from MSigDB (version 7.5.1).}

\item{method}{Name of the method to perform Gene set variation analysis. The
options are: 'gsva', 'ssgea' or 'zscore'. Default value is 'gsva'.}

\item{kcdf}{Character string denoting the kernel to use during the non-parametric
estimation of the cumulative distribution function of expression levels across
samples when method="gsva". By default, "Gaussian" since GEGVIC transforms
raw counts using the vst transformation. Other options are 'Poisson' or 'none'.}

\item{colors}{Character vector indicating the colors of the different groups
to compare. Default values are two: black and orange.}

\item{row.names}{Logical value to determine if row-names are shown in the
heatmap.}

\item{col.names}{Logical value to determine if column-names are shown in the
heatmap.}
}
\value{
Returns a heatmap and the expression values in a form of a matrix.
}
\description{
Performs Gene Set Variation Analysis.
}
\examples{
gsva.res <- ge_single(counts = sample_counts,
                      metadata = sample_metadata,
                      genes_id = 'ensembl_gene_id',
                      response = MSI_status,
                      design = 'MSI_status',
                      biomart = ensembl_biomart_GRCh38_p13,
                      gsva_gmt = 'hallmark',
                      method = 'gsva',
                      kcdf = 'Gaussian',
                      colors = c('orange', 'black'),
                      row.names = TRUE,
                      col.names = TRUE)

}
