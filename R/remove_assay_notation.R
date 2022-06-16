#' Convert DNA sequences with assay notation to regular sequences
#'
#' @param seq DNA sequences in assay notation (i.e. "CCTG\\[T/G\\]TAGG")
#'
#' @return A vector of sequences without assay notation
#' @export
#'
#' @importFrom stringr str_replace
#'
#' @examples
remove_assay_notation <- function(seq) {

  # First make an empty list for the output
  lst <- vector("list", length = length(seq))

  # Loop through the list, getting rid of the extra notation
  for (i in 1:length(seq)) {
    lst[[i]] <- seq[i] %>%
      stringr::str_replace(pattern = "(\\w+)\\[([CATG])\\/[CATG]\\](\\w+)", "\\1\\2\\3")
  }

  # Vectorize the output and return
  return(unlist(lst))
}
