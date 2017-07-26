library(memoise)
ICPerStep <- function(splits, maxIter) ICS(min(splits), max(splits), maxIter)
ICS <- memoise(function(a, b, m) ICSteps(c(rep(1, a), rep(2, b)), maxIter=m))

RandomTrees <- memoise (function(N, n) unclass(rmtree(N, n)))
NamedConstant <- function(X, name) {names(X) <- name; return(X)}

LDFactorial <- memoise (function (x) {
  # memoized version of phangorn::ldfactorial
  x <- (x + 1) / 2
  res <- lgamma(2 * x) - (lgamma(x) + (x - 1) * log(2))
  res
})
LDFact <- memoise(function (x) {
  if (x < 2) return (0) 
  if (x %% 2) {
    LDFactorial(x)
  } else {
    lfactorial(x) - LDFactorial(x - 1L)
  }
})
DFact <- memoise(function (x) exp(LDFact(x)))
DoubleFactorial <- function (x) {
  x[] <- vapply(x, DFact, double(1))
  x
}

NRooted     <- memoise(function (tips, extra=0)  DFact(2 * tips - 3 - extra))
NUnrooted1  <- memoise(function (tips, extra=0)  DFact(2 * tips - 5 - extra))
LnUnrooted1 <- memoise(function (tips, extra=0) LDFact(2 * tips - 5 - extra))
LnRooted    <- memoise(function (tips, extra=0) LDFact(2 * tips - 3 - extra))
N1Spr <- function (n) if (n > 2) 2 * (n - 3) * ((2 * n) - 7) else 0 # Trees exactly one SPR step away. Given by Allen and Steel 2001.

IC1Spr <- function(n) -log2((1+N1Spr(n)) / NUnrooted(n)) # Information content of trees 0 or 1 SPR step from tree with n tips.

LnUnrooted <- function (splits) {
  if ((n.splits <- length(splits)) < 2) return (LnUnrooted1(splits));
  if (n.splits == 2) return (LnRooted(splits[1]) + LnRooted(splits[2]));
  return (LnUnrootedMult(splits))
}
NUnrooted  <- function (splits) {
  if ((n.splits <- length(splits)) < 2) return (NUnrooted1(splits));
  if (n.splits == 2) return (NRooted(splits[1]) *  NRooted(splits[2]))
  return ( NUnrootedMult(splits))
}
LnUnrootedMult <- function (splits) {  # Carter et al. 1990, Theorem 2
  splits <- splits[splits > 0]
  total.tips <- sum(splits)
  LDFact(2 * total.tips - 5) - LDFact(2 * (total.tips - length(splits)) - 1) + sum(vapply(2 * splits - 3, LDFact, double(1)))
}
NUnrootedMult  <- function (splits) {  # Carter et al. 1990, Theorem 2
  splits <- splits[splits > 0]
  total.tips <- sum(splits)
  round(DFact(2 * total.tips - 5) / DFact(2 * (total.tips - length(splits)) - 1) * prod(vapply(2 * splits - 3, DFact, double(1))))
}

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
    tree <- TreeSearch::Postorder(tree)
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

LogisticPoints <- function (x, fitted.model) {
  coefL <- summary(fitted.model)$coef[, 1]
  y <- coefL[1] / (1 + exp((coefL[2] - x) / coefL[3]))
  y
}

Evaluate <- function (tree, data) {
  totalSteps <- TreeSearch::FitchSteps(tree, data)
  chars <- matrix(unlist(data), attr(data, 'nr'))
  ambiguousToken <- which(attr(data, 'allLevels') == "?")
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
  index <- attr(data, 'index')
  total.info <- sum(ic.max[index])
  info.misleading <- sum(info.lost[index])
  proportion.lost <- info.misleading / total.info
  info.needed <- -log(1 / NUnrooted(length(data))) / log(2)
  info.overkill <- total.info / info.needed
  info.retained <- sum(info.used[index])
  signal.noise <- info.retained / info.misleading
  cat("\n", total.info, 'bits, of which', round(info.retained, 2), 'kept,', round(total.info - info.retained, 2), 'lost,', round(info.needed, 2), 'needed.  SNR =', signal.noise, "\n")
  return(c(signal.noise, info.retained/info.needed))
}

#' @returns information content of each extra step, in bits
#' @export
InfoAmounts <- function (data, precision=400000) {
  # The below is simplified from info_extra_step.r::evaluate
  # Assumes no ambiguous tokens & 2 tokens, '1' and '2'
  dataNr <- attr(data, "nr")
  chars <- if (is.null(dim(data))) {
    matrix(c(unlist(data), rep(1, dataNr), rep(2, dataNr)), dataNr)
  } else {
    cbind(data[seq_len(dataNr), ], matrix(c(1, 2), nrow=dataNr, ncol=2, byrow=TRUE))
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
      warning('Concavity function generated by approximation')
    }
    calcIC[seq_along(calculatedP)] <- -log(calculatedP) / log(2)
    calcIC
  }, blankReturn)
  ret[rowSums(ret) > 0, ]
}
