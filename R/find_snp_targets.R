#' Find short microhaplotypes for primer design
#'
#' @param df A data frame with columns 'seq_name' and 'pos'
#' @param span The maximum distance between SNPs when choosing microhaplotypes
#'
#' @return A tibble with columns 'seq_name' and 'target'
#' @export
#'
#' @importFrom dplyr arrange bind_rows
#' @importFrom purrr map
#' @importFrom tibble tibble
#'
#' @examples
find_snp_targets <- function(df, span = 60) {

  lst <- split(df, df$seq_name)

  res_lst <- purrr::map(lst, function(x) {

    seq_id <- x$seq_name[1]
    vec <- x$pos

    if (length(vec) == 1) {
      return(tibble::tibble(seq_name = seq_id, target = paste(vec[1], "1", sep = ",")))
    }

    res <- data.frame(vec = NULL, size = NULL, count = NULL)
    for (i in 1:(length(vec) - 1)) {
      max <- vec[i] + span
      curr_vec <- vec[vec >= vec[i] & vec <= vec[i] + span]
      curr_vec_span <- max(curr_vec) - min(curr_vec) + 1
      res <- rbind(res, data.frame(vec = I(list(curr_vec)), size = curr_vec_span, count = length(curr_vec)))
    }

    res <- dplyr::arrange(res, desc(count), size)
    chosen_vec <- unlist(res$vec[1])
    chosen_range <- res$size[1]

    string <- paste(chosen_vec[1], chosen_range, sep = ",")
    return(tibble::tibble(seq_name = seq_id, target = string))

  })

  final_df <- dplyr::bind_rows(res_lst)

  return(final_df)

}
