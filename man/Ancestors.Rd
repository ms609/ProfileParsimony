\name{Ancestors}
\alias{Ancestors}
\alias{Descendants}
\alias{DoDescendants}
\alias{AllAncestors}
\title{Ancestors and Descendants}
\description{Reports the ancestors and descendants of nodes in a tree.}
\usage{
Ancestors (parent, child, node)
AllAncestors (parent, child)
Descendants (tree, node, ...)
DoDescendants (parent, child, n.tips, node, just.tips = FALSE, just.internal = FALSE,
                           include.ancestor = FALSE)
}
\arguments{
  \item{tree}{An object of class \code{phylo} denoting the tree topology;}
  \item{parent}{the parent node of each edge in a tree, i.e. tree$edge[, 1];}
  \item{child}{the child node of each edge in a tree, i.e. tree$edge[, 2];}
  \item{node}{the number of the node within a tree whose ancestors or descendants should be 
              listed, which might be identified using \code{\link{nodelabels}};}
  \item{n.tips}{the number of tips on the tree;}
  \item{just.tips}{internal nodes are not listed if this is TRUE;}
  \item{just.internal}{tips are not listed if this is TRUE (and just.tips is FALSE);}
  \item{include.ancestor}{includes \code{node} in the output list;}  
  \item{\dots}{other arguments to pass to \code{DoDescendants}.}
}
\details{
  \code{DoDescendants} is more efficient (if less user-friendly) than \code{Descendants}
}
\value{
\itemize{
\item {Ancestors}{A list of all nodes ancestral to \code{node}}
\item {AllAncestors}{A list, each entry comprising a list of all nodes ancestral to each node in turn}
\item {Descendants}{A vector containing descendant nodes in numerical order}
\item {DoDescendants}{A logical vector listing whether each node is a descendant of \code{node}}
}

\author{
Martin R. Smith

Adapted from internal functions in the \pkg{phangorn} package, with performance improvements.
}
\examples{
data(referenceTree)
plot(referenceTree)
nodelabels()
edge   <- referenceTree$edge
parent <- edge[, 1]
child  <- edge[, 2]
Ancestors(parent, child, 36)
AllAncestors(parent, child)[[36]]
Descendants(referenceTree, 36, just.internal=TRUE)
descendant.tips <- DoDescendants(parent, child, length(referenceTree$tip.label), 36, just.tips=TRUE)
referenceTree$tip.label[descendant.tips]

}
\keyword{ tree }