RandomTrees <- memoise (function(N, n) unclass(rmtree(N, n)))
NamedConstant <- function(X, name) {names(X) <- name; return(X)}

ICSteps <- function (char, ambiguous.token = 0, expected.minima = 25, max.iter = 10000) {
  #char <- matrix(c(rep(1, split[1]), rep(2, split[2]), rep(4, split[3]), rep(8, split[4])), ncol=1)
  #data <- phyDat(char, 'USER', 0:3)
  char <- matrix(2 ^ char[char != ambiguous.token], ncol = 1)
  n.char <- length(char)
  rownames(char) <- paste0('t', 1:n.char)
  split <- sort(as.integer(table(char)))
  min.steps <- length(split) - 1
  if (min.steps == 0) return (NamedConstant(1, 0))
  n.no.extra.steps <- n.unrooted.mult(split)
  #n.one.extra.step <- WithOneExtraStep(split)
  proportion.viable <- n.unrooted(n.char) / n.no.extra.steps
  if (proportion.viable == 1) {
    return(NamedConstant(1, min.steps))
  }
  n.iter <- min(max.iter, round(expected.minima * proportion.viable))
  analytic.ic0 <- -log(n.no.extra.steps/n.unrooted(sum(split))) / log(2)
  #analytic.ic1<- -log(n.one.extra.step/n.unrooted(sum(split))) / log(2)
  #analytic.ic1<- -log((n.no.extra.steps + n.one.extra.step)/n.unrooted(sum(split))) / log(2)
  cat(',')
  #cat(c(round(analytic.ic0, 3), 'bits @ 0 extra steps; attempting', n.iter, 'iterations.\n'))
  #cat(c(round(analytic.ic0, 3), 'bits @ 0 extra steps;', round(analytic.ic1, 3), '@ 1; attempting', n.iter, 'iterations.\n'))
  #if (n.iter == max.iter) warning("Truncated number of iterations at max.iter = ", max.iter)
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
  ##table(steps)
  analytic.steps <- n.iter * c(n.no.extra.steps) / n.unrooted(sum(split))
  #analytic.steps <- n.iter * c(n.no.extra.steps, n.one.extra.step) / n.unrooted(sum(split))
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
    backbones <- n.unrooted.mult(backbone.splits)
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
  ic.max <- round(vapply(as.splits, function (split) -log(n.unrooted.mult(split)/n.unrooted(sum(split)))/log(2), double(1)), 12)
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
  info.needed <- -log(1 / n.unrooted(length(data))) / log(2)
  info.overkill <- total.info / info.needed
  info.retained <- sum(info.used[index])
  signal.noise <- info.retained / info.misleading  # With FGH tree: 2.01908.  With pectinate 'basis' tree: 2.033188; calculated with exact values for 1-extra-step: 20.39106.
  cat("\n", total.info, 'bits, of which', round(info.retained, 2), 'kept,', round(total.info - info.retained, 2), 'lost,', round(info.needed, 2), 'needed.  SNR =', signal.noise, "\n")
  return(c(signal.noise, info.retained/info.needed))
}