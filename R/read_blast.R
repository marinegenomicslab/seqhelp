#' Import a BLAST report as a tibble
#'
#' @param file Path to the BLAST report file
#' @param format Currently only supports tab separated reports (outfmt = "6")
#'
#' @return A tibble with BLAST results
#'
#' @note
#' Columns:
#' 1.    qseqid 	 query (e.g., gene) sequence id
#' 2. 	 sseqid 	 subject (e.g., reference genome) sequence id
#' 3. 	 pident 	 percentage of identical matches
#' 4. 	 length 	 alignment length
#' 5. 	 mismatch 	 number of mismatches
#' 6. 	 gapopen 	 number of gap openings
#' 7. 	 qstart 	 start of alignment in query
#' 8. 	 qend 	 end of alignment in query
#' 9. 	 sstart 	 start of alignment in subject
#' 10. 	 send 	 end of alignment in subject
#' 11. 	 evalue 	 expect value
#' 12. 	 bitscore  bit score
#'
#'
#' @export
#'
#' @importFrom attempt stop_if_not
#' @importFrom readr read_tsv
#'
#' @examples
read_blast <- function(file, format = "tab") {

  attempt::stop_if_not(format == "tab", msg = "Need to set format='tab'")

  cl <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

  blast_report <- readr::read_tsv(file,
                    col_names = cl,
                    col_types = readr::cols(
                    qseqid = readr::col_character(),
                    sseqid = readr::col_character(),
                    pident = readr::col_double(),
                    length = readr::col_integer(),
                    mismatch = readr::col_integer(),
                    gapopen = readr::col_integer(),
                    qstart = readr::col_integer(),
                    qend = readr::col_integer(),
                    sstart = readr::col_integer(),
                    send = readr::col_integer(),
                    evalue = readr::col_double(),
                    bitscore = readr::col_double()
                  )
  )

  return(blast_report)

}
