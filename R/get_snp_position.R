#' Locate the position of a SNP in a sequence with assay notation
#'
#' @param seq A vector of sequences in assay notation '\\[C/T\\]'
#'
#' @return An integer representing the location of the SNP in the provided sequences
#' @export
#'
#' @importFrom stringr str_locate
#'
#' @examples
get_snp_position <- function(seq) {

  # First make an empty vector for the output
  lst <- vector("list", length = length(seq))

  # Loop through the sequences, getting the position of the first bracket (the SNP pos)
  for (i in 1:length(seq)) {
    unname(lst[[i]] <- stringr::str_locate(seq[i], "\\[")[,1])
  }

  # Vectorize and return
  return(unname(unlist(lst)))
}
