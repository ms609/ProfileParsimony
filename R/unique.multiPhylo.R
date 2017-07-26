unique.multiPhylo <- function (x, incomparables = FALSE, use.edge.length = FALSE, 
                               use.tip.label = TRUE) {
  n <- length(x)
  if (n < 2) return(x)
  keep <- logical(n)
  keep[1] <- TRUE
  old.index <- seq_len(n)
  for (i in 2:n) {
    already.seen <- FALSE
    comparison.tree <- x[[i]]
    for (j in which(keep)) {
      if (all.equal(x[[j]], comparison.tree, use.edge.length = use.edge.length, 
        use.tip.label = use.tip.label)) {
        already.seen <- TRUE
        old.index[i] <- j
        break
      }
    }
    if (!already.seen) keep[i] <- TRUE
  }
  res <- x[keep]
  attr(res, "old.index") <- old.index
  res
}

