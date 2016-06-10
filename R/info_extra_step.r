library(memoise)
ICPerStep <- function(splits, max.iter) ICS(min(splits), max(splits), max.iter)
ICS <- memoise(function(a, b, m) ICSteps(c(rep(1, a), rep(2, b)), max.iter=m))

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

ICSteps <- function (char, ambiguous.token = 0, expected.minima = 25, max.iter = 10000) {
  char <- matrix(2 ^ char[char != ambiguous.token], ncol = 1)
  n.char <- length(char)
  rownames(char) <- paste0('t', 1:n.char)
  split <- sort(as.integer(table(char)))
  min.steps <- length(split) - 1
  if (min.steps == 0) return (NamedConstant(1, 0))
  n.no.extra.steps <- NUnrootedMult(split)
  #n.one.extra.step <- WithOneExtraStep(split)
  proportion.viable <- NUnrooted(n.char) / n.no.extra.steps
  if (proportion.viable == 1) {
    return(NamedConstant(1, min.steps))
  }
  n.iter <- min(max.iter, round(expected.minima * proportion.viable))
  if (n.iter == max.iter) warning ("Will truncate number of iterations at max.iter = ", max.iter)
  analytic.ic0 <- -log(n.no.extra.steps/NUnrooted(sum(split))) / log(2)
  #analytic.ic1<- -log(n.one.extra.step/NUnrooted(sum(split))) / log(2)
  #analytic.ic1<- -log((n.no.extra.steps + n.one.extra.step)/NUnrooted(sum(split))) / log(2)

  cat(c('  ', signif(analytic.ic0, ceiling(log(max.iter))), 'bits @ 0 extra steps; using', n.iter, 'iterations to estimate cost of further steps.\n'))
  # cat(c(round(analytic.ic0, 3), 'bits @ 0 extra steps;', round(analytic.ic1, 3), '@ 1; attempting', n.iter, 'iterations.\n'))
  trees <- RandomTrees(n.iter, n.char)
  steps <- vapply(trees, function (tree, char) {
    tree <- reorder(tree, "postorder")
    tree.edge <- tree$edge
    tip.label <- tree$tip.label
    parent <- tree.edge[, 1]
    child  <- tree.edge[, 2]
    n.edge <- length(parent)
    max.node <- max(parent)
    stopifnot(max.node == max(tree$edge))
    n.tip <- length(tip.label)
    result <- .Call("FITCH", as.integer(char[tip.label, ]), as.integer(1), 
        as.integer(parent), as.integer(child), as.integer(n.edge), 
        as.double(1), as.integer(max.node), as.integer(n.tip), PACKAGE = "phangorn")
    result
    ret <- result[[1]]
    ret
  }, double(1), char)

  analytic.steps <- n.iter * c(n.no.extra.steps) / NUnrooted(sum(split))
  #analytic.steps <- n.iter * c(n.no.extra.steps, n.one.extra.step) / NUnrooted(sum(split))
  names(analytic.steps) <- min.steps
  #names(analytic.steps) <- min.steps + 0:1
  ##analytic.steps
  tab.steps <- table(steps[steps > (min.steps + 0)])
  #tab.steps <- table(steps[steps > (min.steps + 1)])
  tab.steps <- c(analytic.steps, tab.steps * (n.iter - sum(analytic.steps)) / sum(tab.steps))
  p.steps   <- tab.steps / sum(tab.steps)
  
  return(p.steps)
  #ICSteps  <- -log(cumsum(p.steps)) / log(2)
  #return(rbind(tab.steps, p.steps, ICSteps))
  #summry <- double(length(ICSteps))
  #summry[1:2] <- c(mean(steps), var(steps))
  #return(rbind(tab.steps, p.steps, ICSteps, summry))
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

FitchSite <- function (tree, data) {
  if (class(data) != "phyDat") stop("data must be of class phyDat")
  { #phangorn:::prepareDataFitch
    attrData <- attributes(data)
    lev <- attrData$levels
    l <- length(lev)
    nr <- attrData$nr
    nc <- length(data)
    contrast <- attrData$contrast
    tmp <- contrast %*% 2L^c(0L:(l - 1L))
    tmp <- as.integer(tmp)
    nam <- attrData$names
    attrData$names <- NULL
    data. <- unlist(data, FALSE, FALSE)
    data. <- as.integer(tmp[data.])
    attributes(data.) <- attrData
    attr(data., "dim") <- c(nr, nc)
    dimnames(data.) <- list(NULL, nam)
  }
  { #phangorn:::fit.fitch
    if (is.null(attr(tree, "order")) || attr(tree, "order") == "cladewise") 
        tree <- reorder(tree, "postorder")
    tree.edge <- tree$edge
    node <- tree.edge[, 1]
    edge <- tree.edge[, 2]
    weight <- attrData$weight
    m <- max(tree.edge)
    q <- length(tree$tip)
    result <- .Call("FITCH", data.[, tree$tip.label], as.integer(nr), 
        as.integer(node), as.integer(edge), as.integer(length(edge)), 
        as.double(weight), as.integer(m), as.integer(q), PACKAGE = "phangorn")
    return(result[[2]])
  }
}

Evaluate <- function (tree, data) {
  total.steps <- FitchSite(tree, data)
  chars <- matrix(unlist(data), attr(data, 'nr'))
  ambiguous.token <- which(attr(data, 'allLevels') == "?")
  as.splits <- apply(chars, 1, function (x) {
    ret <- table(x)
    ret[names(ret) != ambiguous.token] 
  })
  if (class(as.splits) == 'matrix') as.splits <- lapply(seq_len(ncol(as.splits)), function(i) as.splits[, i])
  ic.max <- round(vapply(as.splits, function (split) -log(NUnrootedMult(split)/NUnrooted(sum(split)))/log(2), double(1)), 12)
  info.losses <- apply(chars, 1, ICSteps, ambiguous.token=ambiguous.token, max.iter=1000)
  info.amounts <- lapply(info.losses, function(p) {
    #cat(length(p))
    cump <- cumsum(p)
    n.steps <- as.integer(names(p))
    infer <- min(n.steps):max(n.steps)
    infer <- infer[!(infer %in% n.steps)]
    calculated.p <- double(max(n.steps))
    calculated.p[n.steps] <- cump
    if (length(infer)) {
      fitL <- nls(cump ~ SSlogis(n.steps, Asym, xmid, scal))  # Receiving error in p[69
      calculated.p[infer]   <- LogisticPoints(infer, fitL)
    }
    calc.ic <- -log(calculated.p) / log(2)
    calc.ic
  })
  
  info.used <- double(length(total.steps))
  for (i in seq_along(total.steps)) {
    if (total.steps[i] > 0) info.used[i] <- info.amounts[[i]][total.steps[i]]
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

InfoAmounts <- function (data, precision=400000) {
  # The below is simplified from info_extra_step.r::evaluate
  # Assumes no ambiguous tokens & 2 tokens, '1' and '2'
  data.nr <- attr(data, "nr")
  chars <- matrix(c(unlist(data), rep(1, data.nr), rep(2, data.nr)), data.nr) # add 
  splits <- apply(chars, 1, table) - 1
  info.losses <- apply(splits, 2, ICPerStep, max.iter=precision)
  ret <- lapply(info.losses, function(p) {
    cump <- cumsum(p)
    n.steps <- as.integer(names(p))
    infer <- min(n.steps):max(n.steps)
    infer <- infer[!(infer %in% n.steps)]
    calculated.p <- double(max(n.steps))
    calculated.p[n.steps] <- cump
    if (length(infer)) {
      fitL <- nls(cump ~ SSlogis(n.steps, Asym, xmid, scal))
      calculated.p[infer] <- LogisticPoints(infer, fitL)
      warning('Concavity function generated by approximation')
    }
    calc.ic <- -log(calculated.p) / log(2)
    calc.ic
  })
  ret
}
