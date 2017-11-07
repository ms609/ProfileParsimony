#' Profile Parsimony Score
#'
#' Calculate a tree's Profile Parsimony score with a given dataset, after Faith and Trueman (2001)
#'
#' @template treeParam
#' @param dataset Dataset of class \code{profileDat} (see \code{\link{PrepareDataProfile}})
#'                Alternatively, a dataset of class \code{phyDat} can be provided, and will 
#'                be (time-consumingly) converted within the function.
#'
#' @return Zero minus the profile score (because the optimization algorithm assumes that
#'         smaller numbers are better)
#'
#' @references
#'    Faith, D. P. & Trueman, J. W. H. (2001). \cite{Towards an inclusive philosophy for 
#'    phylogenetic inference.} Systematic Biology 50:3, 331-350, doi: 
#'    \href{http://dx.doi.org/10.1080/10635150118627}{10.1080/10635150118627}
#'
#' @examples
#'   data(referenceTree)
#'   data(congreveLamsdellMatrices)
#'   dataset <- PrepareDataFitch(congreveLamsdellMatrices[[42]])
#'   ProfileScore(referenceTree, dataset)
#'
#' @author Martin R. Smith
#'
#' @importFrom TreeSearch Fitch C_Fitch_Steps TipsAreColumns
#' @importFrom TreeSearch TipsAreColumns
#' @keywords tree
#' @export
ProfileScore <- function (tree, dataset) {
  if (class(dataset) == 'phyDat') dataset <- PrepareDataProfile(dataset)
  if (class(dataset) != 'profileDat') {
    stop('Invalid dataset type; prepare dataset with PhyDat() or PrepareDataProfile().')
  }
  at <- attributes(dataset)
  nChar  <- at$nr # strictly, transformation series patterns; these'll be upweighted later
  weight <- at$weight
  steps <- Fitch(tree, dataset, TipData=TipsAreColumns, at, FitchFunction=C_Fitch_Steps)
  info <- at$info.amounts
  nRowInfo <- nrow(info)
  return (-sum(vapply(seq_len(nChar), function (i) {
    stepRow <- max(0L, steps[i] - 1L) + 1L
    return(if (stepRow > nRowInfo) 0 else info[stepRow, i])
  }, double(1)) * weight))
}