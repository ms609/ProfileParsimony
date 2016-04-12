\name{Pratchet}
\alias{Pratchet}
\title{Parsimony ratchet}
\description{This function uses the parsimony ratchet (Nixon 1999) to search for a more parsimonious tree.}
\usage{
Pratchet(tree, data, all = FALSE, outgroup = NULL, pratchiter = 100, searchiter = 5000,
         searchhits = 40, k = 10, track = 0, rearrangements = "NNI", suboptimal = 0, ...) 
}
\arguments{
  \item{tree}{An object of class \code{phyDat} denoting the topology to begin the search from;}
  \item{data}{the output of \code{\link{prepare.data}}, or a matrix in phyDat format
    (which the function will itself pass through \code{prepare.data});}
  \item{all}{Set to TRUE to report all MPTs encountered during the search, perhaps to analyze 
             consensus}
  \item{outgroup}{a vector specifying all tips in the outgroup; if unspecified then identical
                  trees with different roots will be considered unique;}
  \item{pratchiter}{maximum parsimony ratchet iterations to perform;}
  \item{searchiter}{maximum rearrangements to perform on each bootstrap or ratchet iteration;}
  \item{searchhits}{maximum times to hit best score before terminating a tree search within a 
                    pratchet iteration;}
  \item{k}{stop when k ratchet iterations have found the same best score;}
  \item{track}{larger numbers provides more verbose feedback to the user;}
  \item{rearrangements}{method to use when rearranging trees: 'NNI', 'SPR' or 'TBR' !!! 'TBR only' 'SPR only' !!!;}
  \item{suboptimal}{retain trees that are suboptimal by \code{suboptimal} bits;}
  \item{\dots}{other arguments to pass to subsequent functions.}
}
\value{
  This function returns a tree modified by parsimony ratchet iterations,
  retaining the position of the root.
}
\references{
  Nixon, K. C. (1999). \cite{The Parsimony Ratchet, a new method for rapid parsimony analysis.} 
  Cladistics, 15(4), 407-414. doi:\href{http://dx.doi.org/10.1111/j.1096-0031.1999.tb00277.x}{10.1111/j.1096-0031.1999.tb00277.x}
  
  Schliep, K. P. (2011). phangorn: phylogenetic analysis in R. Bioinformatics, 27(4), 592–593.
  doi:\href{http://dx.doi.org/10.1093/bioinformatics/btq706}{10.1093/bioinformatics/btq706}
}

\author{
Martin R. Smith

Adapted from \code{\link{pratchet}} in the \pkg{phangorn} package, which does not preserve the position of the root.
}
\seealso{
\itemize{
\item \code{\link{Pratchet}}
\item \code{\link{TreeSearch}}
}
}
\examples{
\dontrun{
data(congreveLamsdellMatrices)
# Explore trees that are suboptimal by increasing degrees
profile.suboptimal <- c(1e-08, round(0.73^((18:0) - 10), 5))
pp.trees <- lapply(congreveLamsdellMatrices, function (data) 
                   Pratchet(RandomTree(data), PrepareDataFitch(data), 
                            all=TRUE, outgroup='1', suboptimal=max(profile.suboptimal) * 1.1,
                            k=250, searchhits=10, searchiter=50, pratchiter=5000,
                            rearrangements='TBR'))
pp <- vapply(profile.suboptimal, function (so) lapply(pp.trees, function (trees) 
              consensus(trees[Suboptimality(trees) < so])))
plot(pp)
}
}
\keyword{ tree }