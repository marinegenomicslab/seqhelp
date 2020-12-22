#' @keywords internal
"_PACKAGE"

## quiets concerns of R CMD check re: the .'s that appear in pipelines
## https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
if(getRversion() >= "2.15.1") utils::globalVariables(c("."))