Suboptimality <- function (trees, proportional = FALSE) {
  scores <- vapply(trees, attr, double(1), 'score')
  if (proportional) {
    return ((scores - min(scores)) / min(scores))
  } else {
    return(scores - min(scores))
  }
}

Pratchet <- function (tree, data, ParsimonyScorer=ProfileScore, all=FALSE, outgroup=NULL, 
                      pratchiter=100, searchiter=5000, searchhits=40, pratchhits=10, track=0, 
                      rearrangements="NNI", suboptimal=1e-08, ...) {
  if (class(data) == 'phyDat') data <- PrepareDataFitch(data)
  if (class(data) != 'fitchDat') stop("data must be a phyDat object, or the output of PrepareDataFitch(phyDat object).")
  epsilon <- 1e-08
  if (is.null(attr(tree, "score"))) attr(tree, "score") <- ParsimonyScorer(tree, data)
  best.score <- attr(tree, "score")
  if (track >= 0) cat("\n* Initial score:", best.score)
  if (all) {
    null.forest <- vector('list', pratchiter)
    forest <- null.forest
    forest.scores <- rep(NA, pratchiter)
  }

  best.score.hits <- 0
  iterations.completed <- 0
  for (i in 1:pratchiter) {
    if (track >= 0) cat ("\n - Running NNI on bootstrapped dataset. ")
    bstree <- Bootstrap(phy=tree, x=data, maxiter=searchiter, maxhits=searchhits,
                        ParsimonyScorer=ParsimonyScorer,  track=track - 1, ...)
    
    if (track >= 0) cat ("\n - Running", ifelse(is.null(rearrangements), "NNI", rearrangements), "from new candidate tree:")
    if (rearrangements == "TBR") {
      candidate <- TreeSearch(bstree,    data, ParsimonyScorer=ParsimonyScorer, method='TBR', track=track, maxiter=searchiter, maxhits=searchhits, ...)
      candidate <- TreeSearch(candidate, data, ParsimonyScorer=ParsimonyScorer, method='SPR', track=track, maxiter=searchiter, maxhits=searchhits, ...)
      candidate <- TreeSearch(candidate, data, ParsimonyScorer=ParsimonyScorer, method='NNI', track=track, maxiter=searchiter, maxhits=searchhits, ...)
    } else if (rearrangements == "TBR only") {  
      candidate <- TreeSearch(bstree,    data, ParsimonyScorer=ParsimonyScorer, method='TBR', track=track, maxiter=searchiter, maxhits=searchhits, ...)
    } else if (rearrangements == "SPR") {       
      candidate <- TreeSearch(bstree,    data, ParsimonyScorer=ParsimonyScorer, method='SPR', track=track, maxiter=searchiter, maxhits=searchhits, ...)
      candidate <- TreeSearch(candidate, data, ParsimonyScorer=ParsimonyScorer, method='NNI', track=track, maxiter=searchiter, maxhits=searchhits, ...)
    } else if (rearrangements == "SPR only") {  
      candidate <- TreeSearch(bstree,    data, ParsimonyScorer=ParsimonyScorer, method='SPR', track=track, maxiter=searchiter, maxhits=searchhits, ...)
    } else {  
      candidate <- TreeSearch(bstree,    data, ParsimonyScorer=ParsimonyScorer, method='NNI', track=track, maxiter=searchiter, maxhits=searchhits, ...)
    }
    cand.score <- attr(candidate, 'score')
    if ((cand.score + epsilon) < best.score) {
      # New 'best' tree
      if (all) {
        forest[[i]] <- if (is.null(outgroup)) candidate else Root(candidate, outgroup)
        forest.scores[i] <- cand.score
      }
      tree <- candidate
      best.score <- cand.score
      best.score.hits <- 1
    } else if (best.score + epsilon > cand.score) { # i.e. best == cand, allowing for floating point error
      best.score.hits <- best.score.hits + 1
      tree <- candidate
      if (all) {
        forest[[i]] <- if (is.null(outgroup)) candidate else Root(candidate, outgroup)
        forest.scores[i] <- cand.score
      }
    } else if (cand.score < (best.score + suboptimal) && all) {
      forest[[i]] <- if (is.null(outgroup)) candidate else Root(candidate, outgroup)
      forest.scores[i] <- cand.score
    }
    if (track >= 0) cat("\n* Best score after", i, "/", pratchiter, "pratchet iterations:", best.score, "( hit", best.score.hits, "/", pratchhits, ")")
    if (best.score.hits >= pratchhits) {
      iterations.completed <- i
      break()
    }
  } # end for
  if (iterations.completed == 0) iterations.completed <- pratchiter
  if (track >= 0) cat ("\nCompleted parsimony ratchet after", iterations.completed, "iterations with score", best.score, "\n")
   
  if (all) {
    keepers <- !is.na(forest.scores) & forest.scores < best.score + suboptimal
    forest.scores <- forest.scores[keepers]
    forest <- forest[keepers]
    if (length(forest) > 1) {
      class(forest) <- 'multiPhylo'
      ret <- unique(forest)
    } else if (length(forest) == 1) {
      class(forest) <- 'phylo'
      ret <- forest
    } else {
      stop('No trees!? Is suboptimal set to a sensible (positive) value?')
    }
    scores.unique <- vapply(ret, attr, double(1), 'score')
    cat('Found', sum(scores.unique == min(scores.unique)), 'unique MPTs and', length(ret) - sum(scores.unique == min(scores.unique)), 'suboptimal trees.\n')
    if (is.null(outgroup)) warning('"outgroup" not specified, so some "unique" trees may have same topology but distinct roots.')
  } else {
    ret <- tree
    attr(ret, 'hits') <- NULL
  }
  return (ret)
}
  
Bootstrap <- function (phy, x, maxiter, maxhits, ParsimonyScorer = ProfileScore, track=1, ...) {
## Simplified version of phangorn::bootstrap.phyDat, with bs=1 and multicore=FALSE
  at <- attributes(x)
  weight <- at$weight
  v <- rep(1:length(weight), weight)
  BS <- tabulate(sample(v, replace=TRUE), length(weight)) 
  keep <- BS > 0
  ind <- which(keep)
  x <- x[ind, ]
  attr(x, 'weight') <- BS[ind]
  attr(x, 'min.steps') <- at$min.steps[keep]
  attr(x, 'info.amounts') <- at$info.amounts[keep, ]
  attr(x, 'unique.tokens') <- at$unique.tokens[keep]
  attr(x, 'sa.weights') <- at$sa.weights[keep]
  attr(x, 'nr') <- length(ind)
  attr(x, 'inapp.level') <- at$inapp.level
  attr(phy, 'score') <- NULL
  class(x) <- 'fitchDat'
  res <- TreeSearch(phy, x, ParsimonyScorer=ParsimonyScorer, method='NNI', maxiter=maxiter,
                    maxhits=maxhits, track=track-1, ...)
  attr(res, 'score') <- NULL
  attr(res, 'hits') <- NULL
  res
}

TreeSearch <- function (tree, data, ParsimonyScorer = ProfileScore,  method = 'NNI', 
                        maxiter = 100, maxhits = 20, forest.size = 1,
                        cluster = NULL, track = 1, ...) {
  tree$edge.length <- NULL # Edge lengths are not supported
  attr(tree, 'hits') <- 1
  if (exists("forest.size") && length(forest.size) && forest.size > 1) {
    forest <- empty.forest <- vector('list', forest.size)
    forest[[1]] <- tree
  } else {
    forest.size <- 1 
  }
  if (is.null(attr(tree, 'score'))) attr(tree, 'score') <- ParsimonyScorer(tree, data)
  best.score <- attr(tree, 'score')
  if (track > 0) cat("\n  - Performing", method, "search.  Initial score:", best.score)
  Rearrange <- switch(method, 'TBR' = TBR, 'SPR' = SPR, 'NNI' = NNI)
  return.single <- !(forest.size > 1)
  
  for (iter in 1:maxiter) {
    trees <- RearrangeTree(tree, data, Rearrange, ParsimonyScorer, min.score=best.score,
                           return.single=return.single, iter=iter, cluster=cluster, track=track)
    iter.score <- attr(trees, 'score')
    if (length(forest.size) && forest.size > 1) {
      hits <- attr(trees, 'hits')
      if (iter.score == best.score) {
        forest[(hits-length(trees)+1L):hits] <- trees
        tree <- sample(forest[1:hits], 1)[[1]]
        attr(tree, 'score') <- iter.score
        attr(tree, 'hits') <- hits
      } else if (iter.score < best.score) {
        best.score <- iter.score
        forest <- empty.forest
        forest[1:hits] <- trees
        tree <- sample(trees, 1)[[1]]
        attr(tree, 'score') <- iter.score
        attr(tree, 'hits') <- hits
      }      
    } else {
      if (iter.score <= best.score) {
        best.score <- iter.score
        tree <- trees
      }
    }
    if (attr(trees, 'hits') >= maxhits) break
  }
  if (track > 0) cat("\n  - Final score", attr(tree, 'score'), "found", attr(tree, 'hits'), "times after", iter, method, "iterations\n")  
  if (forest.size > 1) {
    if (hits < forest.size) forest <- forest[-((hits+1):forest.size)]
    attr(forest, 'hits') <- hits
    attr(forest, 'score') <- best.score
    return(unique(forest))
  } else {
    return(tree)
  }
}

# Adapted from phangorn:::nni; speed increased
NNI <- function (tree) {
  n      <- sample(tree$Nnode - 1L, 1L)
  edge   <- tree$edge
  parent <- edge[, 1L]
  child  <- edge[, 2L]
  root   <- min(parent)
  n.tips <- root - 1L
  ind    <- which(child > n.tips)[n]
  p1     <- parent[ind]
  p2     <- child[ind]
  ind1   <- which(parent == p1)
  ind1   <- ind1[ind1 != ind][1L]
  ind2   <- which(parent == p2)[sample(2L,1L)]
  # Perform the switch
  child[c(ind1, ind2)] <- child[c(ind2, ind1)]
    
  # Reorder pruning - modified from phangorn:::reorderPruning  
  max.node <- max(parent)
  n.edges  <- length(parent)
  neworder <- .C("C_reorder", as.integer(parent), as.integer(child), as.integer(n.edges), 
                 as.integer(max.node), integer(n.edges), as.integer(root-1L), PACKAGE = "phangorn")[[5]]    
  tree$edge <- matrix(c(parent, child), ncol=2)[neworder,]
  tree$edge.length <- tree$edge.length[neworder]
  attr(tree, "order") <- "pruningwise"
  Renumber(tree)
}

SPR <- function(tree) {
  tip.label <- tree$tip.label
  nTips <- length(tip.label)
  edge  <- tree$edge; parent <- edge[,1L]; child <- edge[,2L]
  nEdge <- length(child)
  root  <- nTips + 1L
  if (nTips < 4L) stop ('must be >3 tips for SPR rearrangement!')
  pruning.candidates <- seq(nEdge + 1L)[-root]
  repeat {
    prune.node <- sample(pruning.candidates, 1)
    moving.subnodes <- c(prune.node, which(DoDescendants(parent, child, nTips, prune.node)))
    moving.nodes <- c(prune.parent <- parent[child==prune.node], moving.subnodes)
    dont.graft.here <- c(moving.nodes, child[parent==prune.parent])
    graft.node <- c(pruning.candidates[!pruning.candidates %in% dont.graft.here])
    if (length(graft.node) > 1) graft.node <- sample(graft.node, 1)
    if (any(graft.node)) break;
    pruning.candidates <- pruning.candidates[-match(prune.node, pruning.candidates)]
    if (!any(pruning.candidates)) stop('No place to graft pruned tree')
  } 
  
  graft.edge   <- match(graft.node, child)
  graft.parent <- parent[graft.edge]
  graft.child  <-  child[graft.edge]
  prune.edge   <- match(prune.node, child)
  parent.duplicate <- parent
  parent.duplicate[prune.edge] <- NA
  sister.edge  <- match(prune.parent, parent.duplicate)
  if (prune.parent == root) {
    new.root <- child[parent==root]
    new.root <- new.root[new.root != prune.node]
    edge[sister.edge, 2L] <- edge[graft.edge, 2L]
    edge[graft.edge, 2L] <- root
    new.root.spots <- edge==new.root
    edge[edge == root] <- new.root
    edge[new.root.spots] <- root
  } else {
    leading.edge <- match(prune.parent, child)
    edge[c(leading.edge, sister.edge, graft.edge), 2] <- edge[c(sister.edge, graft.edge, leading.edge), 2]
  }
  
  reordered.edge <- .C('order_edges', as.integer(edge[,1]), as.integer(edge[,2]),
                       as.integer(nTips-1L), as.integer(nEdge), PACKAGE='ProfileParsimony')
  numbered.edge <- .C('number_nodes', as.integer(reordered.edge[[1]]), 
                      as.integer(reordered.edge[[2]]), as.integer(root), as.integer(nEdge), 
                      PACKAGE='ProfileParsimony')
  tree$edge <- matrix(c(numbered.edge[[1]], numbered.edge[[2]]), ncol=2)
  tree
}

TBR <- function(tree, edge.to.break=NULL) {
  nTips <- tree$Nnode + 1
  if (nTips < 3) return (tree)
  tree.edge <- tree$edge
  tree.parent <- tree.edge[,1]
  tree.child <- tree.edge[,2]
  if (nTips == 3) return (Root(tree, sample(tree.child[tree.parent==max(tree.parent)], 1L)))
  all.nodes <- 1:(2*(nTips-1))
  root <- nTips + 1
  if (is.null(edge.to.break)) edge.to.break <- sample(2L:nrow(tree.edge), 1L) # Only include one root edge
  subtree.root <- tree.child[edge.to.break]
  stump <- if (subtree.root <= nTips) {
    DropTipNoSubtree(tree, subtree.root, root.edge=1)
  } else {
    in.crown <- DoDescendants(tree.parent, tree.child, nTips, subtree.root, just.tips=TRUE)
    DropTipNoSubtree(tree, which(in.crown), root.edge=1)
  }
  stump.len <- dim(stump$edge)[1]
  crown <- ExtractClade(tree, subtree.root)  # ~ 2x faster than DropTip
  crown.edge <- crown$edge
  crown.len <- dim(crown.edge)[1]  
  if (crown.len > 1) {
    if (crown.len == 2) {
      rerooted.crown <- crown
    } else {
      crown.parent <- crown.edge[,1]
      crown.child <- crown.edge[,2]
      crown.nNode <- crown$Nnode
      crown.tips <- crown.nNode + 1L
      crown.root <- min(crown.parent)
      new.root.candidates <- crown.child[-1] # Include existing root once only
      new.root.node <- sample(new.root.candidates, 1L)
      new.outgroup <- if (new.root.node <= crown.tips) {
        new.root.node 
      } else {
        which(DoDescendants(crown.parent, crown.child, crown.tips,
                            new.root.node, just.tips=TRUE))
      }
      rerooted.crown <- Root(crown, new.outgroup)
    }
    rerooted.crown$root.edge <- 1L
    if (stump.len > 1) {
      bind.location <- sample(seq_len(stump.len), 1L)
      ret <- BindTree(stump, rerooted.crown, position=1, where=bind.location)
    } else {
      ret <- AddTip(rerooted.crown, crown.root, stump$tip.label)
    }
  } else {
    bind.location <- stump$edge[sample(2L:stump.len, 1L), 2L]
    ret <- AddTip(stump, bind.location, crown$tip.label)
  }
  Renumber(ret)
}

Fitch <- function (tree, data, at = NULL) {
  if (is.null(at)) at <- attributes(data)
  n.char  <- at$nr # strictly, transformation series patterns; these'll be upweighted later
  weight <- at$weight
  info <- at$info.amounts
  if (is.null(at$order) || at$order == "cladewise") tree <- reorder(tree, "postorder")
  tree.edge <- tree$edge
  parent <- tree.edge[, 1]
  child <- tree.edge[, 2]
  tip.label <- tree$tip.label
  n.edge <- length(parent)
  max.node <- parent[1] #max(parent)
  n.tip <- length(tip.label)
  n.node <- max.node - n.tip
  inapp <- at$inapp.level
  parent.of <- parent[match((n.tip + 2L):max.node, child )]
  allNodes <- (n.tip + 1L):max.node
  child.of <- child [c(match(allNodes, parent),
                       length(parent) + 1L - match(allNodes, rev(parent)))]
  fitch <- .Call("FITCH", data[, tree$tip.label], as.integer(n.char),
        as.integer(parent), as.integer(child), as.integer(n.edge),
        as.double(weight), as.integer(max.node), as.integer(n.tip), PACKAGE='phangorn')
  return(fitch[[2]])
#
#  Future support for inapplicable data to be added here:
#  
#  nLevel <- length(at$level)
#  powers.of.2 <- 2L^c(0L:(nLevel - 1L))
#  inapp.level <- which(at$levels == "-")
#  applicable.tokens <- setdiff(powers.of.2, 2^(inapp.level - 1))
#  fitch <- .Call("FITCHINAPP", data[, tip.label], as.integer(n.char), as.integer(parent),
#                 as.integer(child), as.integer(parent.of), as.integer(child.of),
#                 as.integer(n.edge), as.integer(n.node), as.double(weight),
#                 as.integer(max.node), as.integer(n.tip), as.integer(inapp),
#                 PACKAGE='inapplicable') 
}

ProfileScore <- function (tree, data) {
    # Data
  if (class(data) == 'phyDat') data <- PrepareDataFitch(data)
  if (class(data) != 'fitchDat') stop('Invalid data type; try ProfileScore(tree, data <- 
                                      PrepareDataFitch(valid.phyDat.object)).')
  at <- attributes(data)
  n.char  <- at$nr # strictly, transformation series patterns; these'll be upweighted later
  weight <- at$weight
  steps <- Fitch(tree, data, at)
  # Return a negative rather than positive value because algorithms assume that 
  # smaller numbers are better
  return(-sum(info[(steps - 1) * n.char + seq_len(n.char)] * weight))
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

SuccessiveApproximations <- function (tree, data, outgroup = NULL, k = 3, max.succiter = 20,
                                      pratchhits = 100, searchhits = 50, searchiter = 500,
                                      pratchiter = 5000, track = 0, suboptimal = 0.1) {
  
  if (k < 1) stop ('k should be at least 1, see Farris 1969 p.379')
  attr(data, 'sa.weights') <- rep(1, length(attr(data, 'weight')))
  collect.suboptimal <- suboptimal > 0
  
  max.node <- max(tree$edge[, 1])
  n.tip <- length(tree$tip.label)
  n.node <- max.node - n.tip
  bests <- vector('list', max.succiter + 1)
  bests.consensus <- vector('list', max.succiter + 1)
  best <- bests[[1]] <- bests.consensus[[1]] <- Root(tree, outgroup)
  for (i in seq_len(max.succiter) + 1) {
    if (track > 0) cat('\nSuccessive Approximations Iteration', i - 1)
    attr(best, 'score') <- NULL
    if (suboptimal > 0) {
      suboptimal.search <- suboptimal * sum(attr(data, 'sa.weights') * attr(data, 'weight'))
    }
    trees <- Pratchet(best, data, ParsimonyScorer = SuccessiveWeights, all = collect.suboptimal, 
                           suboptimal=suboptimal.search,    rearrangements='NNI',
                           pratchhits=pratchhits, searchhits=searchhits, searchiter=searchiter, 
                           pratchiter=pratchiter, outgroup = outgroup, track=track - 1)
    trees <- unique(trees)
    bests[[i]] <- trees
    suboptimality <- Suboptimality(trees)
    bests.consensus[[i]] <- consensus(trees[suboptimality == 0])
    if (all.equal(bests.consensus[[i]], bests.consensus[[i - 1]])) return(bests[2:i])
    best <- trees[suboptimality == 0][[1]]
    l.i <- Fitch(best, data)
    p.i <- l.i / (n.node - 1)
    w.i <- ((p.i)^-k) - 1
    attr(data, 'sa.weights') <- w.i
  }
  cat('Stability not reached.')
  return(bests)
}
