Ancestors <- function (parent, child, node) {
  if (length(node) == 1) {
    pvector <- numeric(max(parent))
    pvector[child] <- parent
    anc <- function(pvector, node) {
      res <- numeric(0)
      repeat {
        anc <- pvector[node]
        if (anc == 0) 
            break
        res <- c(res, anc)
        node <- anc
      }
      res
    }
    return(anc(pvector, node))
  } else AllAncestors(parent, child)[node]
}

AllAncestors <- function (parent, child) {
  res <- vector("list", max(parent))
  for (i in seq_along(parent)) {
    pa <- parent[i]
    res[[child[i]]] <- c(pa, res[[pa]])
  }
  res
}

Descendants <- function (tree, node, ...) {
# ARGUMENTS:
#   "tree", object of class phylo
#   "node", number of an internal node
# RETURN:
#   vector containing descendant nodes in numerical order
  nTip <- length(tree$tip.label)
  edge <- tree$edge
  edge1 <- edge[,1]
  edge2 <- edge[,2]
  return (which(DoDescendants(edge1, edge2, nTip, node, ...)))
}

DoDescendants <- function (parent, child, n.tips, node, just.tips = FALSE, just.internal = FALSE,
                           include.ancestor = FALSE) {
  # ARGUMENTS:
  #   "parent", parent nodes: from tree$edge[,1]
  #   "child", parent nodes: from tree$edge[,2]
  #   "node", number of an internal node
  #   "just.tips", should return value include all nodes or just tips?
  # RETURN:
  #   vector containing descendant nodes in numerical order
  is.descendant <- blank <- logical((n.tips * 2) - 1)
  if (include.ancestor) is.descendant[node] <- TRUE;
  node.children <- function (node, is.descendant) {
    nc <- child[parent %in% node]
    is.descendant[nc] <- TRUE
    if (length(nc)) is.descendant <- node.children(nc, is.descendant)
    is.descendant
  }
  is.descendant <- node.children(node, is.descendant)
  if (just.tips) return (is.descendant[1:n.tips]) else if (just.internal) is.descendant[1:n.tips] <- FALSE 
  return (is.descendant)
}
