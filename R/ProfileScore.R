#' Profile Score
#'
#' @template treeParam
#' @param data Dataset of class \code{fitchDat} 
#' @returns Zero minus the profile score (because the optimization algorithm assumes that
#' smaller numbers are better)
#' @importFrom TreeSearch FitchSteps
#' @importFrom TreeSearch TipsAreColumns
#' @keywords internal
#' @export
ProfileScore <- function (tree, data) {
  if (class(data) == 'phyDat') data <- PrepareDataFitch(data)
  if (class(data) != 'fitchDat') {
    stop('Invalid data type; prepare data with PhyDat() or PrepareDataFitch().')
  }
  at <- attributes(data)
  nChar  <- at$nr # strictly, transformation series patterns; these'll be upweighted later
  weight <- at$weight
  steps <- TreeSearch::FitchSteps(tree, data, TreeSearch::TipsAreColumns, at)
  info <- at$info.amounts
  nRowInfo <- nrow(info)
  return (-sum(vapply(seq_len(nChar), function (i) {
    stepRow <- max(0L, steps[i] - 1L) + 1L
    return(if (stepRow > nRowInfo) 0 else info[stepRow, i])
  }, double(1)) * weight))
}