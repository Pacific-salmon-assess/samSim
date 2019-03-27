#' ungroupRowwiseDF
#'
#' Simple helper function to undo rowwise grouping.
#'
#' @param data Dataframe that is grouped rowwise.
#' @return Returns a dataframe that is no longer grouped rowwise.
#'
#' @export
ungroupRowwiseDF <- function(x) {
  class(x) <- c( "tbl_df", "data.frame")
  x
}
