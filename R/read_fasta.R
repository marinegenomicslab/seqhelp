#' Import a FASTA file into a tibble
#'
#' @param fasta Path to the input FASTA file
#'
#' @return A tibble with columns 'seq_id' and 'seq'
#' @export
#'
#' @importFrom seqinr getName getSequence read.fasta
#' @importFrom tibble tibble
#'
#' @examples
read_fasta <- function(fasta) {

  # Read the FASTA file with seqinr and create a tibble
  input_fa <- seqinr::read.fasta(fasta)
  seq <- toupper(unlist(seqinr::getSequence(input_fa, as.string = TRUE)))
  seq_ids <- seqinr::getName(input_fa)
  fa_tbl <- tibble::tibble(seq_id = seq_ids, seq = seq)

  return(fa_tbl)

}
