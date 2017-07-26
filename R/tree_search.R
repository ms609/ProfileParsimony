Suboptimality <- function (trees, proportional = FALSE) {
  scores <- vapply(trees, attr, double(1), 'score')
  if (proportional) {
    return ((scores - min(scores)) / min(scores))
  } else {
    return(scores - min(scores))
  }
}

#' Profile Score
#'
#' @template treeParam
#' @param data Dataset of class \code{fitchDat} 
#' @returns Zero minus the profile score (because the optimization algorithm assumes that
#' smaller numbers are better)
#' @keywords internal
#' @export
ProfileScore <- function (tree, data) {
  if (class(data) == 'phyDat') data <- PrepareDataFitch(data)
  if (class(data) != 'fitchDat') {
    stop('Invalid data type; prepare data with PhyDat() or PrepareDataFitch().')
  }
  at <- attributes(data)
  n.char  <- at$nr # strictly, transformation series patterns; these'll be upweighted later
  weight <- at$weight
  steps <- Fitch(tree, data, at)
  info <- at$info.amounts
  return(-sum(info[max(0, (steps - 1)) * n.char + seq_len(n.char)] * weight))
}

SuccessiveWeights <- function(tree, data) {
  # Data
  if (class(data) == 'phyDat') data <- PrepareDataFitch(data)
  if (class(data) != 'fitchDat') {
    stop('Invalid data type; prepare data with PhyDat() or PrepareDataFitch().')
  }
  at <- attributes(data)
  weight <- at$weight
  sa.weights <- at$sa.weights
  if (is.null(sa.weights)) sa.weights <- rep(1, length(weight))
  steps <- Fitch(tree, data, at)
  return(sum(steps * sa.weights * weight))
}
