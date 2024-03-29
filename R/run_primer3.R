#' Design primers for a sequence
#'
#' @param df A dataframe with at least the following columns: seq, seq_name, target_string
#' seq: A character string containing the sequence
#' seq_name: The name of the sequence
#' target_string: String to identify the target region - takes the form "target,length"
#' @param amplicon_length Vector of length three with (min, opt, max) of amplicon product size
#' @param melting_temp Vector of length three with (min, opt, max) of melting temperature (Tm)
#' @param primer_length Vector of length three with (min, opt, max) of primer length
#' @param settings_file Optional path to a file of predetermined settings
#' @param settings Optional named list of additional settings
#'
#' @return
#' @export
#'
#' @importFrom dplyr everything filter mutate pull select row_number
#' @importFrom glue glue
#' @importFrom purrr walk2
#' @importFrom readr read_delim write_lines
#' @importFrom stringr str_replace_all
#' @importFrom tidyr extract pivot_wider
#'
#' @examples
run_primer3 <- function(df,
                        amplicon_length = c(110, 150, 190),
                        melting_temp = c(60, 60, 65),
                        primer_length = c(20, 20, 25),
                        settings_file = NULL,
                        settings = NULL) {


  attempt::stop_if_not(all(c("seq", "seq_name", "target") %in% colnames(df)),
                       msg = "Input dataframe must have the following columns: seq, seq_name, target")


  #browser()

  # Find the executable paths
  src <- file.path(system.file(package = "seqhelp"), "exec", "primer3_core")
  therm_params <- file.path(system.file(package = "seqhelp"), "primer3_config")

  # Optionally use a custom settings file
  if (!is.null(settings_file)) {
    ps_settings <- settings_file
  } else { # Otherwise, use the default
    primer3_settings <- file.path(system.file(package = "seqhelp"), "primer3_settings.txt")
  }

  # Loop through the rows

  results <- purrr::map(1:nrow(df), function(i) {

    seq <- df$seq[i]
    seq_name <- df$seq_name[i]
    target_string <- df$target[i]

    # Create tempfiles
    primer3_input <- tempfile()
    primer3_output <- tempfile()

    # Write the input file
    readr::write_lines(sprintf("SEQUENCE_ID=%s", seq_name), primer3_input)
    readr::write_lines(sprintf("SEQUENCE_TEMPLATE=%s", as.character(seq)), primer3_input, append = TRUE)
    readr::write_lines("PRIMER_TASK=generic", primer3_input, append = TRUE)
    readr::write_lines("PRIMER_PICK_LEFT_PRIMER=1", primer3_input, append = TRUE)
    readr::write_lines("PRIMER_PICK_INTERNAL_OLIGO=0", primer3_input, append = TRUE)
    readr::write_lines("PRIMER_PICK_RIGHT_PRIMER=1", primer3_input, append = TRUE)
    readr::write_lines("PRIMER_EXPLAIN_FLAG=1", primer3_input, append = TRUE)
    readr::write_lines("PRIMER_PAIR_MAX_DIFF_TM=3", primer3_input, append = TRUE)
    readr::write_lines(sprintf("PRIMER_MIN_TM=%s", melting_temp[1]), primer3_input, append = TRUE)
    readr::write_lines(sprintf("PRIMER_OPT_TM=%s", melting_temp[2]), primer3_input, append = TRUE)
    readr::write_lines(sprintf("PRIMER_MAX_TM=%s", melting_temp[3]), primer3_input, append = TRUE)
    readr::write_lines(sprintf("PRIMER_MIN_SIZE=%s", primer_length[1]), primer3_input, append = TRUE)
    readr::write_lines(sprintf("PRIMER_OPT_SIZE=%s", primer_length[2]), primer3_input, append = TRUE)
    readr::write_lines(sprintf("PRIMER_MAX_SIZE=%s", primer_length[3]), primer3_input, append = TRUE)
    readr::write_lines(sprintf("PRIMER_PRODUCT_SIZE_RANGE=%s-%s", amplicon_length[1], amplicon_length[3]), primer3_input, append = TRUE)
    readr::write_lines(sprintf("PRIMER_PRODUCT_OPT_SIZE=%s", amplicon_length[2]), primer3_input, append = TRUE)
    readr::write_lines(sprintf("PRIMER_THERMODYNAMIC_PARAMETERS_PATH=%s", therm_params), primer3_input, append = TRUE)

    if (! is.null(target_string)) {
      readr::write_lines(sprintf("SEQUENCE_TARGET=%s", target_string), primer3_input, append = TRUE)
    }

    if (! is.null(settings)) {
        purrr::walk2(settings, names(.), ., function(x, y) {
          readr::write_lines(sprintf("%s=%s", x, y), primer3_input, append = TRUE)
      })
    }

    readr::write_lines("=", primer3_input, append = TRUE)

    # Run Primer3
    try(system(glue::glue("{src} {primer3_input} -p3_settings_file {primer3_settings} > {primer3_output}")))

    # Read in the output
    primer_res <- readr::read_delim(primer3_output, delim = '=', col_names = c("param", "value"), col_types = readr::cols())

    if ("PRIMER_ERROR" %in% primer_res$param) {

      return(list(summary = NULL, results = NULL))
    }


    # Extract the last summary row number
    summary_row_end <- primer_res %>%
      dplyr::mutate(row = dplyr::row_number()) %>%
      dplyr::filter(param == "PRIMER_PAIR_NUM_RETURNED") %>%
      dplyr::pull(row)



    if(length(summary_row_end) == 0) {browser()}

    # Create a summary data frame from the primer candidates
    summary_tbl_raw <- primer_res[1:summary_row_end,]

    summary_tbl <- tibble::tibble(locus = summary_tbl_raw$value[1],
                                  summary_left = summary_tbl_raw$value[19],
                                  summary_right = summary_tbl_raw$value[20],
                                  summary_pair = summary_tbl_raw$value[21],
                                  pairs_returned = summary_tbl_raw$value[25])



    if (primer_res$value[primer_res$param == "PRIMER_PAIR_NUM_RETURNED"][1] == 0) {
      return(list(summary = summary_tbl, results = NULL))
    }

    primer_tbl_long <- primer_res[(summary_row_end + 1):nrow(primer_res),]

    # Convert the primer data frame to a format with one row per primer pair and clean it up
    primer_tbl_wide <- tidyr::extract(primer_tbl_long, param, "pair_idx", regex = "_(\\d+)", remove = FALSE)
    primer_tbl_wide <- dplyr::filter(primer_tbl_wide, !is.na(pair_idx))
    primer_tbl_wide <- dplyr::mutate(primer_tbl_wide, col = stringr::str_replace_all(param, "_\\d+_", "_"))
    primer_tbl_wide <- dplyr::mutate(primer_tbl_wide, col = stringr::str_replace_all(col, "_\\d+", ""))
    primer_tbl_wide <- tidyr::pivot_wider(primer_tbl_wide, names_from = col, values_from = value, id_cols = pair_idx)
    primer_tbl_wide <- tidyr::extract(primer_tbl_wide, PRIMER_LEFT, c("PRIMER_LEFT_START", "PRIMER_LEFT_LENGTH"), "(\\d+),(\\d+)")
    primer_tbl_wide <- tidyr::extract(primer_tbl_wide, PRIMER_RIGHT, c("PRIMER_RIGHT_START", "PRIMER_RIGHT_LENGTH"), "(\\d+),(\\d+)")
    primer_tbl_wide <- dplyr::mutate(primer_tbl_wide, seq_name = seq_name)
    primer_tbl_wide <- dplyr::select(primer_tbl_wide, seq_name, rank = pair_idx, dplyr::everything())

    return(list(summary = summary_tbl, results = primer_tbl_wide))

  })

  cand_seqs <- dplyr::bind_rows(purrr::transpose(results)$results)

  summaries <- dplyr::bind_rows(purrr::transpose(results)$summary)

  # Clean up the summary tibble

  summaries <- dplyr::mutate(summaries, status = dplyr::if_else(.data$pairs_returned > 0, "PASSED", "FAILED"))

  # TODO: More informative summaries of failed primers
  summaries <- tidyr::extract(summaries,
                              col = c("summary_left"),
                              into = c("left_considered", "left_fails", "left_passed"),
                              regex = "considered (\\d+), (.*)ok (\\d+)")
  summaries <- tidyr::extract(summaries,
                              col = c("summary_right"),
                              into = c("right_considered", "right_fails", "right_passed"),
                              regex = "considered (\\d+), (.*)ok (\\d+)")
  summaries <- tidyr::extract(summaries,
                              col = c("summary_pair"),
                              into = c("pairs_considered", "pairs_fails", "pairs_passed"),
                              regex = "considered (\\d+), (.*)ok (\\d+)")


  return(list(primers = cand_seqs, summary = summaries))


}

