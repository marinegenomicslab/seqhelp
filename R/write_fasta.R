#' Write a FASTA file
#'
#' @param file The path of the FASTA file to create
#' @param names Vector of sequence names (to become the FASTA ID)
#' @param seq Vector of sequences
#'
#' @return Nothing
#' @export
#'
#' @importFrom seqinr write.fasta
#'
#' @examples
write_fasta <- function(file, names, seq) {

  # Write a FASTA file with seqinr
  seqinr::write.fasta(sequences = as.list(seq),
                      names = names,
                      file.out = file,
                      as.string = TRUE)

}
