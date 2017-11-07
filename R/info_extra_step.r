#' Information content per step
#' @export
ICPerStep <- function(splits, maxIter) ICS(min(splits), max(splits), maxIter)

#' Information content per step
#' @importFrom R.cache addMemoization
#' @export
ICS <- addMemoization(function(a, b, m) ICSteps(c(rep(1, a), rep(2, b)), maxIter=m))

#' Random trees
#' @importFrom memoise memoise
#' @export
RandomTrees <- memoise (function(N, n) unclass(rmtree(N, n)))
NamedConstant <- function(X, name) {names(X) <- name; return(X)}

#' Log double factorial
#' 
#' Memoised version of phangorn's \code{\link[phangorn]{ldfactorial}}
#' @importFrom phangorn ldfactorial
#' @importFrom memoise memoise
#' @export
LDFactorial <- memoise(ldfactorial)

#' Log double factorial
#' Handles odd and even inout numbers
#' @param x a positive integer
#' @importFrom memoise memoise
#' @export
LDFact <- memoise(function (x) {
  if (x < 2) return (0) 
  if (x %% 2) {
    LDFactorial(x)
  } else {
    lfactorial(x) - LDFactorial(x - 1L)
  }
})

#' Double factorial
#' @param x a positive integer
#' @export
DFact <- memoise(function (x) exp(LDFact(x)))

#' @describeIn DFact Accepts a vector as input
#' @param ints a vector of integers
#' @export
DoubleFactorial <- function (ints) {
  ints[] <- vapply(ints, DFact, double(1))
  ints
}

#' Number of rooted/unrooted trees
#' These functions return the number of rooted or unrooted trees consistent with a given pattern
#'  of splits.
#'
#' 
#' Functions starting N return the number of rooted or unrooted trees, functions starting Ln
#' provide the log of this number.  Calculations follow Carter et al. 1990, Theorem 2.
#'
#' @param tips The number of tips.
#' @param extra the number of points at which another branch cannot be added.
#' @param splits vector listing the number of taxa in each tree bipartition.
#'
#' @author Martin R. Smith
#' 
#' @references 
#' Carter, M., Hendy, M., & Penny, D. (1990). \cite{On the distribution of lengths of evolutionary 
#' trees.} SIAM Journal on Discrete Mathematics, 3(1), 38-47. doi:
#' \href{http://doi.org/10.1137/0403005}{10.1137/0403005}
#' @examples
#'   NRooted(10)
#'   NUnrooted(10)
#'   LnRooted(10)
#'   LnUnrooted(10)
#'   # Number of trees consistent with a character whose states are 00000 11111 222
#'   NUnrootedMult(c(5,5,3))
#' 
#' @export
NRooted     <- memoise(function (tips, extra=0)  DFact(2 * tips - 3 - extra))
#' @describeIn NRooted Number of unrooted trees
#' @export
NUnrooted1  <- memoise(function (tips, extra=0)  DFact(2 * tips - 5 - extra))
#' @describeIn NRooted  Log Number of unrooted trees
#' @export
LnUnrooted1 <- memoise(function (tips, extra=0) LDFact(2 * tips - 5 - extra))
#' @describeIn NRooted  Log Number of rooted trees
#' @export
LnRooted    <- memoise(function (tips, extra=0) LDFact(2 * tips - 3 - extra))
#' Number of trees on SPR step away
#' Formula given by Given by Allen and Steel 2001.
#' @references ALLEN, B. L. and STEEL, M. 2001. Subtree transfer operations and their 
#'             induced metrics on evolutionary trees. \emph{Annals of Combinatorics},
#'             5, 1--15. <doi:10.1007/s00026-001-8006-8>
#' @export
N1Spr <- function (n) if (n > 2) 2 * (n - 3) * ((2 * n) - 7) else 0 

#' @describeIn N1Spr Information content of trees 0 or 1 SPR step from tree with n tips.
#' @export
IC1Spr <- function(n) -log2((1+N1Spr(n)) / NUnrooted(n)) 

#' @describeIn NRooted Log number of unrooted trees
#' @template splitsParam
#' @export
LnUnrooted <- function (splits) {
  if ((n.splits <- length(splits)) < 2) return (LnUnrooted1(splits));
  if (n.splits == 2) return (LnRooted(splits[1]) + LnRooted(splits[2]));
  return (LnUnrootedMult(splits))
}
#' @describeIn NRooted Number of unrooted trees
#' @export
NUnrooted  <- function (splits) {
  if ((n.splits <- length(splits)) < 2) return (NUnrooted1(splits));
  if (n.splits == 2) return (NRooted(splits[1]) *  NRooted(splits[2]))
  return ( NUnrootedMult(splits))
}
#' @describeIn NRooted Log unrooted mult
#' @references CARTER, M., HENDY, M., PENNY, D., SZEKELY, L. A. and WORMALD, N. C. 1990.
#'             On the distribution of lengths of evolutionary trees. 
#'             \emph{SIAM Journal on Discrete Mathematics}, 3, 38--47.
#' @export
LnUnrootedMult <- function (splits) {  # Carter et al. 1990, Theorem 2
  splits <- splits[splits > 0]
  total.tips <- sum(splits)
  LDFact(2 * total.tips - 5) - LDFact(2 * (total.tips - length(splits)) - 1) + sum(vapply(2 * splits - 3, LDFact, double(1)))
}
#' @describeIn NRooted Number of unrooted trees (mult)
#' @export
NUnrootedMult  <- function (splits) {  # Carter et al. 1990, Theorem 2
  splits <- splits[splits > 0]
  total.tips <- sum(splits)
  round(DFact(2 * total.tips - 5) / DFact(2 * (total.tips - length(splits)) - 1) * prod(vapply(2 * splits - 3, DFact, double(1))))
}

#' Information Content Steps
#'
#'   This function estimates the information content of a character \code{char} when e extra steps
#'   are present, for all possible values of e.
#' 
#' @param char The character in question.
#' @param ambiguous.token Which token, if any, corresponds to the ambigious token (?)
#'                          (not yet fully implemented).
#' @param expected.minima sample enough trees that the rarest step counts is expected to be seen
#'                          at least this many times.
#' @param max.iter Maximum iterations to conduct.
#' 
#' Calculates the number of trees consistent with the character having \emph{e} extra steps, where
#' e ranges from its minimum possible value (i.e. number of different tokens minus one) to its
#' maximum.  The number of trees with no extra steps can be calculated exactly; the number of trees
#' with more additional steps must be approximated.
#' The function samples \code{n.iter} trees, or enough trees that the trees with the minimum number 
#' of steps will be recovered at least \code{expected.minima} times, in order to obtain precise 
#' results.
#' 
#' @author{
#' Martin R. Smith
#' }
#' @references{
#' 
#' Faith, D. P. & Trueman, J. W. H. (2001). \cite{Towards an inclusive philosophy for phylogenetic
#' inference.} Systematic Biology 50:3, 331-350, doi:
#' \href{http://dx.doi.org/10.1080/10635150118627}{10.1080/10635150118627}
#' 
#' }
#' @keywords tree
#' 
#' @examples{
#'   # A character that is present in ten taxa and absent in five
#'   character <- c(rep(1, 10), rep(2, 5))
#'   ICSteps (character)
#' }
#' @importFrom TreeSearch Postorder C_Fitch_Score
#' @export
ICSteps <- function (char, ambiguousToken = 0, expectedMinima = 25, maxIter = 10000) {
  char <- matrix(2 ^ char[char != ambiguousToken], ncol = 1)
  rownames(char) <- paste0('t', seq_along(char))
  charLen <- length(char)
  
  split <- sort(as.integer(table(char)))
  minSteps <- length(split) - 1
  if (minSteps == 0) return (NamedConstant(1, 0))
  
  nNoExtraSteps <- NUnrootedMult(split)
  #nOneExtraStep <- WithOneExtraStep(split)
  proportionViable <- NUnrooted(charLen) / nNoExtraSteps
  if (proportionViable == 1) return(NamedConstant(1, minSteps))
  
  nIter <- min(maxIter, round(expectedMinima * proportionViable))
  if (nIter == maxIter) warning ("Will truncate number of iterations at maxIter = ", maxIter)
  analyticIc0 <- -log(nNoExtraSteps/NUnrooted(sum(split))) / log(2)
  #analyticIc1<- -log(nOneExtraStep/NUnrooted(sum(split))) / log(2)
  #analyticIc1<- -log((nNoExtraSteps + nOneExtraStep)/NUnrooted(sum(split))) / log(2)

  cat('  Token count', split, "=", signif(analyticIc0, ceiling(log10(maxIter))), 'bits @ 0 extra steps; simulating', nIter, 'trees to estimate cost of further steps.\n')
  # cat(c(round(analyticIc0, 3), 'bits @ 0 extra steps;', round(analyticIc1, 3), '@ 1; attempting', nIter, 'iterations.\n'))
  trees <- RandomTrees(nIter, charLen)  ## TODO make more efficient by randomising trees that are already in postorder
  
  nEdge   <- 2L * charLen - 2L
  maxNode <- 2L * charLen - 1L
  
    
  steps <- vapply(trees, function (tree, char) {
    tree <- Postorder(tree)
    treeEdge <- tree$edge
    tipLabel <- tree$tip.label
    parent <- treeEdge[, 1]
    child  <- treeEdge[, 2]
    return(C_Fitch_Score(char[tipLabel, ], nChar=1L, parent, child, nEdge, weight=1L,
        maxNode, nTip=charLen))    
  }, double(1), char)

  analyticSteps <- nIter * c(nNoExtraSteps) / NUnrooted(sum(split))
  #analyticSteps <- nIter * c(nNoExtraSteps, nOneExtraStep) / NUnrooted(sum(split))
  names(analyticSteps) <- minSteps
  #names(analyticSteps) <- minSteps + 0:1
  ##analyticSteps
  tabSteps <- table(steps[steps > (minSteps + 0)]) # Quicker than table(steps)[-0]
  #tabSteps <- table(steps[steps > (minSteps + 1)])
  tabSteps <- c(analyticSteps, tabSteps * (nIter - sum(analyticSteps)) / sum(tabSteps))
  pSteps <- tabSteps / sum(tabSteps)
  
  return(pSteps)
  #ICSteps  <- -log(cumsum(pSteps)) / log(2)
  #return(rbind(tabSteps, pSteps, ICSteps))
  #summry <- double(length(ICSteps))
  #summry[1:2] <- c(mean(steps), var(steps))
  #return(rbind(tabSteps, pSteps, ICSteps, summry))
}

#' Numer of trees with one extra step
#' @export
WithOneExtraStep <- function (split) {
  # Ignore singletons, which can be added at the end...
  split.with.splittables <- split[split > 1]
  if (length(split.with.splittables) < 2) return (0)
  
  # TODO test split: 1 1 2 4, split: 2 2 4
  sum(vapply(seq_along(split), function (omit) {
    backbone.splits <- split[-omit]
    omitted.tips <- split[omit]
    if (omitted.tips < 2) return (0)
    backbone.tips <- sum(backbone.splits)
    backbones <- NUnrootedMult(backbone.splits)
    backbone.edges <- max(0, 2 * backbone.tips - 3)
    backbone.attachments <- backbone.edges * (backbone.edges - 1)
    prod(sum( # Branch unambiguously split along first group
      vapply(1:(omitted.tips - 1), function (first.group) { # For each way of splitting up the omitted tips, e.g. 1|16, 2|15, 3|14, etc
        choose(omitted.tips, first.group) * 
        n.rooted(first.group) * n.rooted(omitted.tips - first.group)
      }, double(1))
    ) / 2, backbone.attachments, backbones) + prod(
    # Second group added adjacent to first group, thus new edge could belong to the backbone or the omitted tip group
    sum(vapply(1:(omitted.tips - 1), function (first.group) { # For each way of splitting up the omitted tips, e.g. 1|16, 2|15, 3|14, etc
        choose(omitted.tips, first.group) * 
        n.rooted(first.group) * n.rooted(omitted.tips - first.group) # backbone tips have already been split - when we selected a branch
      }, double(1))) / 2,
    backbones,
    backbone.edges,
    2 # left or right of group addition location
    / 2 # Will be counted again when 'added group' becomes the 'backbone group'
    )
    
  }, double(1))
  )
}

#' @export
LogisticPoints <- function (x, fitted.model) {
  coefL <- summary(fitted.model)$coef[, 1]
  y <- coefL[1] / (1 + exp((coefL[2] - x) / coefL[3]))
  y
}

#' Evaluate tree
#' @template treeParam
#' @template datasetParam
#' @importFrom TreeSearch Fitch C_Fitch_Steps
#' @export
Evaluate <- function (tree, dataset) {
  totalSteps <- Fitch(tree, dataset, FitchFunction = C_Fitch_Steps)
  chars <- matrix(unlist(dataset), attr(dataset, 'nr'))
  ambiguousToken <- which(attr(dataset, 'allLevels') == "?")
  as.splits <- apply(chars, 1, function (x) {
    ret <- table(x)
    ret[names(ret) != ambiguousToken] 
  })
  if (class(as.splits) == 'matrix') as.splits <- lapply(seq_len(ncol(as.splits)), function(i) as.splits[, i])
  ic.max <- round(vapply(as.splits, function (split) -log(NUnrootedMult(split)/NUnrooted(sum(split)))/log(2), double(1)), 12)
  infoLosses <- apply(chars, 1, ICSteps, ambiguousToken=ambiguousToken, maxIter=1000)
  infoAmounts <- lapply(infoLosses, function(p) {
    #cat(length(p))
    cumP <- cumsum(p)
    nSteps <- as.integer(names(p))
    infer <- min(nSteps):max(nSteps)
    infer <- infer[!(infer %in% nSteps)]
    calculatedP <- double(max(nSteps))
    calculatedP[nSteps] <- cumP
    if (length(infer)) {
      fitL <- nls(cumP ~ SSlogis(nSteps, Asym, xmid, scal))  # Receiving error in p[69
      calculatedP[infer]   <- LogisticPoints(infer, fitL)
    }
    calcIC <- -log(calculatedP) / log(2)
    calcIC
  })
  
  info.used <- double(length(totalSteps))
  for (i in seq_along(totalSteps)) {
    if (totalSteps[i] > 0) info.used[i] <- infoAmounts[[i]][totalSteps[i]]
  }
  info.lost <- round(ic.max - info.used, 13)
  index <- attr(dataset, 'index')
  total.info <- sum(ic.max[index])
  info.misleading <- sum(info.lost[index])
  proportion.lost <- info.misleading / total.info
  info.needed <- -log(1 / NUnrooted(length(dataset))) / log(2)
  info.overkill <- total.info / info.needed
  info.retained <- sum(info.used[index])
  signal.noise <- info.retained / info.misleading
  cat("\n", total.info, 'bits, of which', round(info.retained, 2), 'kept,', round(total.info - info.retained, 2), 'lost,', round(info.needed, 2), 'needed.  SNR =', signal.noise, "\n")
  # Return:
  c(signal.noise, info.retained/info.needed)
}

#' Amount of information in each character
#'
#' @template datasetParam
#' @param precision number of random trees to generate when calculating Profile curves
#' @template warnParam
#'
#' @return information content of each extra step, in bits
#' 
#' @author Martin R. Smith
#'
#' @export
InfoAmounts <- function (dataset, precision=4e+05, warn=TRUE) {
  # The below is simplified from info_extra_step.r::evaluate
  # Assumes no ambiguous tokens & 2 tokens, '1' and '2'
  dataNr <- attr(dataset, "nr")
  chars <- if (is.null(dim(dataset))) {
    matrix(c(unlist(dataset), rep(1, dataNr), rep(2, dataNr)), dataNr)
  } else {
    cbind(dataset[seq_len(dataNr), ], matrix(c(1, 2), nrow=dataNr, ncol=2, byrow=TRUE))
  }
  splits <- apply(chars, 1, table) - 1
  infoLosses <- apply(splits, 2, ICPerStep, maxIter=precision)
  
  blankReturn <- double(max(colSums(splits)))
  ret <- vapply(infoLosses, function(p) {
    calcIC <- blankReturn
    cumP <- cumsum(p)
    nSteps <- as.integer(names(p))
    infer <- min(nSteps):max(nSteps)
    infer <- infer[!(infer %in% nSteps)]
    calculatedP <- double(max(nSteps))
    calculatedP[nSteps] <- cumP
    if (length(infer)) {
      fitL <- nls(cumP ~ SSlogis(nSteps, Asym, xmid, scal))
      calculatedP[infer] <- LogisticPoints(infer, fitL)
      if (warn) warning('Concavity function generated by approximation')
    }
    calcIC[seq_along(calculatedP)] <- -log(calculatedP) / log(2)
    calcIC
  }, blankReturn)
  ret[rowSums(ret) > 0, ]
}
