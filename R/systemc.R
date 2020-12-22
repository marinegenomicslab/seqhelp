#' Run a system call inside a conda environment
#'
#' @param command The system command (Linux only)
#' @param env The conda environment
#' @param ignore_stdout Ignore STDOUT (TRUE/FALSE)
#' @param ignore_stderr Ignore STDERR (TRUE/FALSE)
#' @param conda_sh_path The path to the conda.sh
#'
#' @return The return value from the system call
#' @export
#'
#' @importFrom glue glue
#'
#' @examples
systemc <- function(command, env = CONDA_ENV, ignore_stdout = FALSE, ignore_stderr = FALSE, conda_sh_path = CONDA_PATH) {
  com <- glue::glue("bash -c '. {conda_sh_path}; conda activate {env}; {command}'")
  system(com, ignore.stdout = ignore_stdout, ignore.stderr = ignore_stderr)
}
