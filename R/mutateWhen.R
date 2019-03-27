#' mutate_when
#'
#' Simple helper function to allow conditional application of mutate.
#'
#' @param data Dataframe that meets standard dplyr requirements.
#' @return Returns a dataframe subsetted as requested.
#'
#' @examples
#' mtcars %>% mutate_when(
#'    mpg > 22,    list(cyl = 100),
#'    disp == 160, list(cyl = 200)
#' )
#' @export
mutate_when <- function(data, ...) {
  dots <- eval(substitute(alist(...)))
  for (i in seq(1, length(dots), by = 2)) {
    condition <- eval(dots[[i]], envir = data)
    mutations <- eval(dots[[i + 1]], envir = data[condition, , drop = FALSE])
    data[condition, names(mutations)] <- mutations
  }
  data
}
