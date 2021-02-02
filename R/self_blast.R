#' BLAST-compare a set of sequences to identify regions of similarity
#'
#' @param seq A vector of sequences
#' @param id A vector of sequence names
#' @param conda_env The conda environment
#' @param conda_sh_path The path to the conda.sh
#'
#' @return A tibble with regions of similarity
#' @export
#'
#' @importFrom glue glue
#' @importFrom dplyr arrange desc filter
#' @importFrom rlang .data
#'
#' @examples
self_blast <- function(seq, id, conda_env = CONDA_ENV, conda_sh_path = CONDA_PATH) {

  # Make some temporary files to hold the FASTA file and BLAST DB
  tmp_dir <- tempdir()
  tmp_fasta <- glue::glue(tmp_dir, "tmp.fasta", .sep = "/")
  tmp_db <- glue::glue(tmp_dir, "tmp", .sep = "/")
  tmp_out <- glue::glue(tmp_dir, "tmp_out.blast", .sep = "/")

  # Write the FASTA file
  write_fasta(file = tmp_fasta, seq = seq, names = id)

  # Set up the BLAST DB and run the self-BLAST
  invisible({
    systemc(glue::glue("makeblastdb -in {tmp_fasta} -dbtype nucl -out {tmp_db}"), env = conda_env, ignore_stdout = TRUE, ignore_stderr = TRUE, conda_sh_path = conda_sh_path)
    systemc(glue::glue("blastn -query {tmp_fasta} -db {tmp_db} -outfmt 6 -task blastn-short > {tmp_out}"), env = conda_env, conda_sh_path = conda_sh_path)
  })

  # Read in the results
  blast_results <- read_blast(tmp_out)

  # Remove self-hits and sort by bitscore
  blast_results <- dplyr::filter(blast_results, .data$qseqid != .data$sseqid)
  blast_results <- dplyr::arrange(blast_results, dplyr::desc(.data$bitscore))

  # Return the BLAST results
  return(blast_results)

}
