#' @title ic_raw_to_tpm
#'
#' @description Transforms RNA-seq raw counts to TPM (Transcript Per kilobase
#' Million).
#'
#' @param counts Data frame that contains gene expression data as raw counts.
#' @param genes_id Name of the column that contains gene identifiers. Should be
#' one of the following:'entrezgene_id', 'ensembl_gene_id' or 'hgnc_symbol'.
#' @param biomart Data frame containing a biomaRt query with the following
#' attributes: ensembl_gene_id, hgnc_symbol, entrezgene_id, transcript_length,
#' refseq_mrna.
#'
#' @return Returns a matrix containing expression counts as TPM with HGNC gene
#' symbols as rownames and samples identifiers as colnames.
#'
#' @export
#'
#' @import dplyr
#' @import tibble
#'
#' @examples
#' tpm <- ic_raw_to_tpm(counts = input_ge_module,
#'                      genes_id = 'entrezgene_id',
#'                      biomart = ensembl_biomart_GRCh38_p13)

ic_raw_to_tpm <- function(counts,
                          genes_id,
                          biomart) {

    # Get the length (in Kb) for each gene
    # Filter transcripts that have refseq_mrna ID in BiomaRt object
    biomart_length <- biomart %>%
        # filter genes that have refseq Id and gene symbol
        dplyr::filter(refseq_mrna != '' & hgnc_symbol != '') %>%
        # group_by gene symbol
        dplyr::group_by(hgnc_symbol) %>%
        # Calculate median length in Kb for each gene
        dplyr::summarise(median_kb_transcript_length = median(transcript_length) / 1000)

    # If data is annotated with entrezgene_id convert them to character
    if (genes_id == 'entrezgene_id'){
        raw.count <- counts %>%
            dplyr::mutate(entrezgene_id = as.character(entrezgene_id))
    } else {
        raw.count <- counts
    }

    # Process the raw.counts
    raw.count <- raw.count %>%
            # Join the data frame with the biomaRt table
        dplyr::inner_join(x = .,
                          y = biomart %>%
                              mutate(entrezgene_id = as.character(entrezgene_id)),
                          by = genes_id) %>%
        # From all annotation columns keep only hgnc symbol column
        dplyr::select(hgnc_symbol,
                      transcript_length,
                      everything(),
                      -c(entrezgene_id, ensembl_gene_id,
                         transcript_length, refseq_mrna)) %>%
        # Filter those missing gene symbols
        dplyr::filter(hgnc_symbol != '') %>%
        # Remove duplicated genes
        dplyr::distinct(hgnc_symbol, .keep_all = TRUE) %>%
        # Get the rows that are present in the biomart_length object
        dplyr::semi_join(., biomart_length, by = 'hgnc_symbol')

    # Get also for biomart_length object those genes in raw.counts so we can
    ## can join both objects
    biomart_length <- biomart_length %>%
        dplyr::semi_join(., raw.count, by = 'hgnc_symbol')

    # Calculate RPK (reads per kilobase)
    rpk <- raw.count %>%
        # Reorder the rows as in biomart_length
        dplyr::arrange(match(x = hgnc_symbol,
                      table = biomart_length$hgnc_symbol)) %>%
        tibble::column_to_rownames('hgnc_symbol')


    # Create a function to calculate TPM from RPK
    tpm_calc <- function(rpk,len) {
        x <- rpk/len
        return(t(t(x)*1e6/colSums(x)))
    }

    # Calculate TPM from RPK
    tpm <- tpm_calc(rpk = rpk,
                    len = biomart_length$median_kb_transcript_length)

    return(tpm)
}
