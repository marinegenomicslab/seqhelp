#' Fetch a nucleotide sequence from NCBI or a FASTA file
#'
#' @param accession Sequence accession number from NCBI nuccore
#' @param fasta Path to the fasta file
#' @param seq_name The name of the sequence in the FASTA file
#' @param seq_start Start coordinate
#' @param seq_end End coordinate
#'
#' @return The requested sequence (character string)
#' @export
#'
#' @importFrom glue glue
#' @importFrom IRanges IRanges IRangesList
#' @importFrom Rsamtools FaFile scanFa
#'
#'
#' @examples
fetch_seq <- function(accession = NULL, fasta = NULL, seq_name = NULL, seq_start, seq_end) {

  attempt::stop_if(!is.null(accession) & !is.null(fasta), msg = "Need to specify either an accession or FASTA file, not both!")

  if (!is.null(fasta)) {
    irange <- IRanges::IRanges(start = c(seq_start), end = c(seq_end), names = c(seq_name))
    irange_list <- IRanges::IRangesList(irange)
    names(irange_list) <- seq_name
    faidx <- Rsamtools::FaFile(fasta)
    grange <- Rsamtools::scanFa(faidx, param = irange_list)
    seq <- as.character(grange$seq)
  }


  if (!is.null(accession)) {

    # System call
    cmd <-  glue::glue("wget -q -O - 'https://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={acc}&seq_start={start}&seq_stop={end}&rettype=fasta'")

    seq <- system(cmd, intern = TRUE)
    seq <- seq[2:(length(seq) - 1)]
    seq <- paste0(seq, collapse = "")

    # Check to see if the sequence was returned properly
    if (nchar(seq) != as.numeric(seq_end) - as.numeric(seq_start) + 1) {
      cat(nchar(seq))
      cat(as.numeric(seq_end) - as.numeric(seq_start) + 1)
      stop("Error fetching sequence")
    }

  }

  return(seq)
}
