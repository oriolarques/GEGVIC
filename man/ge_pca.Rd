% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ge_pca.R
\name{ge_pca}
\alias{ge_pca}
\title{ge_pca}
\usage{
ge_pca(
  counts,
  genes_id,
  metadata,
  response,
  design,
  colors = c("orange", "black")
)
}
\arguments{
\item{counts}{Data frame that contains gene expression data as raw counts.}

\item{genes_id}{Name of the column that contains gene identifiers. Should be
one of the following:'entrezgene_id', 'ensembl_gene_id' or 'hgnc_symbol'.}

\item{metadata}{Data frame that contains supporting variables to the data.}

\item{response}{Unquoted name of the variable indicating the groups to analyse.}

\item{design}{Variables in the design formula in the form of: 'Var1 + Var2 + ... Var_n'.}

\item{colors}{Character vector indicating the colors of the different groups
to compare. Default values are two: black and orange.}
}
\value{
Returns a ggplot object.
}
\description{
Creates a PCA plot from gene expression data using the VST
transformation from the plotPCA function from the DESeq2 package.
}
\examples{
ge_pca(counts = sample_counts,
       genes_id = 'ensembl_gene_id',
       metadata = sample_metadata,
       response = MSI_status,
       design = 'MSI_status',
       colors = c('orange', 'black'))

}
