\name{RearrangeTree}
\alias{RearrangeTree}
\title{Rearrange phylogenetic tree}
\description{Performs one tree rearrangement of the specified type.}
\usage{
RearrangeTree(tree, data, Rearrange = NNI, ParsimonyScorer = ProfileScore, min.score = NULL, 
              return.single = TRUE, iter = "<unknown>", cluster = NULL, trace = 0)
}
\arguments{
  \item{tree}{a rooted bifurcating phylogenetic tree with the desired outgroup, and the attributes
    \code{pscore}, the tree's parsimony score, and 
    \code{hits}, the number of times the best score has been hit in the calling function;
  }
  \item{data}{a data matrix defined by \code{\link{prepare.data}};}
  \item{Rearrange}{a rearrangement function: probably one of 
    \code{\link{NNI}}, \code{\link{SPR}} or \code{\link{TBR}};}
  \item{ParsimonyScorer}{a function to generate optimality score;}
  \item{min.score}{trees longer than \code{min.score}, probably the score of the starting tree,
    will be discarded;}
  \item{return.single}{returns all trees if \kbd{FALSE} or a randomly selected tree if \kbd{TRUE};}
  \item{iter}{iteration number of calling function, for reporting to user via stdout only;}
  \item{cluster}{a cluster, to accelerate searches on multicore machines (not currently implemented);}
  \item{trace}{determines how much information to output to screen.}
}
\value{This function performs one iteration of \code{Rearrange} on the \code{tree} (on each core
 in the cluster, if used).  If a cluster is used, the most parsimonious tree is returned - 
 or most parsimonious trees, if multiple are found and \code{return.single = FALSE}.
 Each tree is annotated with attributes 'hits' (the number of time that tree was recovered) and 
 'score' (parsimony score of the tree under \code{ParsimonyScorer}).
 as described for argument \code{tree}.}
\author{Martin R. Smith}
\seealso{
  \itemize{
    \item \code{\link{NNI}}
    \item \code{\link{SPR}}
    \item \code{\link{TBR}}
  }
}
\examples{
data('referenceTree')
data('congreveLamsdellMatrices')
dat <- congreveLamsdellMatrices[[52]]
parsimony(referenceTree, dat)
tree <- RearrangeTree(referenceTree, dat, SPR, parsimony)
par(mfrow=c(1, 2))
plot(referenceTree)
plot(tree)
attributes(tree)
}