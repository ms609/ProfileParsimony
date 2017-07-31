#' Profile Score
#'
#' Calculate a tree's Profile Parsimony score, after Faith and Trueman (2001)
#'
#' @template treeParam
#' @param data Dataset of class \code{phyDat} or (preferably) \code{profileDat} 
#'             (see \code{\link{PrepareDataProfile}})
#' @return Zero minus the profile score (because the optimization algorithm assumes that
#'         smaller numbers are better)
#' @importFrom TreeSearch FitchSteps
#' @importFrom TreeSearch TipsAreColumns
#'
#' @author Martin R. Smith
#'
#' @export
ProfileScore <- function (tree, data) {
  if (class(data) == 'phyDat') data <- PrepareDataProfile(data)
  if (class(data) != 'fitchDat') {
    stop('Invalid data type; prepare data with PhyDat() or PrepareDataProfile().')
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