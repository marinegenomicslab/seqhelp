#' Extend assay sequences using a reference genome
#'
#' @param tbl A table of assays with columns: locus, acc, gstart, gend, reverse, pos, ref, alt, orig_seq
#' @param len Number of bases to extend on either side of the original sequence
#'
#' @return A tibble with the extended sequence
#' @export
#'
#' @importFrom attempt stop_if_not
#' @importFrom Biostrings DNAString reverseComplement
#' @importFrom dplyr bind_rows
#' @importFrom glue glue
#' @importFrom purrr map
#' @importFrom stringr str_length str_locate str_sub
#' @importFrom tibble tibble
#'
#' @examples
extend_seq_range <- function(tbl, len) {

  attempt::stop_if_not(all(c("locus", "acc", "gstart", "gend", "reverse", "pos", "ref", "alt", "orig_seq") %in% colnames(tbl)),
                       msg = "Input table must have the following columns: locus, acc, gstart, gend, reverse, pos, ref, alt, orig_seq")

  tbl %>%
    split(tbl$locus) %>%
    purrr::map(function(x) {

      tag_length <- 8
      ext_seq <- get_genome_seq(acc = x$acc, start = x$gstart - len, end = x$gend + len) %>%
        Biostrings::DNAString()

      if (isTRUE(x$reverse)) {
        ext_seq <- Biostrings::reverseComplement(ext_seq)
        new_start <- x$gend + len
        new_end <- x$gstart - len
      } else {
        new_end <- x$gend + len
        new_start <- x$gstart - len
      }

      # Get a 17 bp tag for the SNP
      tag <- stringr::str_sub(x$orig_seq, start = x$pos - tag_length, end = x$pos + tag_length)
      tag_mat <- Biostrings::matchPattern(tag, ext_seq, max.mismatch = 1)

      if (length(tag_mat@ranges) == 1) {

        # Get the beginning and end of the final string
        string_begin <- stringr::str_sub(as.character(ext_seq), 1, tag_mat@ranges@start - 1)
        string_end <- stringr::str_sub(as.character(ext_seq), tag_mat@ranges@start + tag_mat@ranges@width, stringr::str_length(as.character(ext_seq)))

        # Build the SNP string
        snp_string_start <- stringr::str_sub(as.character(ext_seq), tag_mat@ranges@start, tag_mat@ranges@start + tag_length - 1)
        snp_site <- glue::glue("[{x$ref}/{x$alt}]")
        snp_string_end <- stringr::str_sub(as.character(ext_seq), tag_mat@ranges@start + tag_length + 1, tag_mat@ranges@start + (tag_length*2))
        snp_string <- paste0(snp_string_start, snp_site, snp_string_end)

        # Put all of the pieces together
        final_seq <- paste0(string_begin, snp_string, string_end)

        # Get the new position of the SNP
        new_pos <- stringr::str_locate(final_seq, "\\[")[,1] %>%
          unname()

        return(tibble::tibble(locus = x$locus,
                      acc = x$acc,
                      gstart = new_start,
                      gend = new_end,
                      reverse = x$reverse,
                      pos = new_pos,
                      ref = x$ref,
                      alt = x$alt,
                      ext_seq = final_seq)
        )


      } else {
        #browser()
        return(tibble::tibble(locus = x$locus,
                      acc = x$acc,
                      gstart = new_start,
                      gend = new_end,
                      reverse = x$reverse,
                      pos = NA,
                      ref = x$ref,
                      alt = x$alt,
                      ext_seq = NA)
        )
      }

    }) %>%
    dplyr::bind_rows()

}
