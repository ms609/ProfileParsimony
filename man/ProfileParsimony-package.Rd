\name{ProfileParsimony-package}
\alias{ProfileParsimony-package}
\alias{ProfileParsimony}
\docType{package}
\title{Phylogenetic Inference using Profile Parsimony}
\description{
This package implements the Profile Parsimony approach of Faith & Trueman (2001), which finds 
the tree that is most faithful to the information contained within a given dataset.
It also provides functions to rearrange trees whilst preserving the outgroup. 
The package also contains heuristic search methods to locate the most parsimonious tree.
}
\details{
This package calculates the parsimony score on phylogenetic trees.  It can also be used to find
the most parsimonious tree.  The \emph{Examples} section below provides a step-by-step guide to 
using this package on your own data.

Character data are read from a \emph{restricted} \acronym{NEXUS} format using the R function
\code{\link{read.nexus.data}}; see the latter function's \link[=read.nexus.data]{documentation} for details.

\acronym{NEXUS} files can be edited in any standard text editor,
 for example \href{http://notepad-plus-plus.org/}{Notepad++}.
 
A notable annoyance is that the parser cannot interpret curly braces, for example {01}.  
Ambiguous tokens of this nature should be replaced with a separate character, for example 
'A' to denote {01}, 'B' to denote {12}, perhaps using a search-and-replace operation in your 
favourite text editor.

The algorithm is currently implemented only for binary characters.  This shortcoming will 
be addressed soon.

The package includes 100 simulated datasets (from Congreve & Lamsdell 2016), used in Smith (201X)
as a test of the efficacy of Profile Parsimony.
}

\author{
Martin R. Smith
}
\references{
Congreve, C. R., & Lamsdell, J. C. (2016). \cite{Implied weighting and its utility in 
  palaeontological datasets: a study using modelled phylogenetic matrices.} Palaeontology, 
  online ahead of print. doi:10.1111/pala.12236

Congreve, C. R., & Lamsdell, J. C. (2016). Data from: Implied weighting and its utility in 
palaeontological datasets: a study using modelled phylogenetic matrices. Dryad Digital
Repository, \href{http://dx.doi.org/10.5061/dryad.7dq0j}{doi:10.5061/dryad.7dq0j}

Faith, D. P. & Trueman, J. W. H. (2001). \cite{Towards an inclusive philosophy for phylogenetic
inference.} Systematic Biology 50:3, 331-350, doi:
\href{http://dx.doi.org/10.1080/10635150118627}{10.1080/10635150118627}

Smith, M. R. (201X). \cite{Trade-Offs between Resolution and Accuracy under Different Methods of
Phylogenetic Inference.} Syst. Biol. 
}
\keyword{ package }
\keyword{ tree }
\seealso{
\itemize{
\item \code{\link[ape:ape-package]{ape}}
\item \code{\link[phangorn:phangorn-package]{phangorn}}
}}
\examples{
\dontrun{

## Walkthrough of package functions

## To use this script on your own data, launch R, and type (or copy-paste) the following text
## into the R console.  Lines starting with '#' are comments and do not need to be copied.


## To install the package for the first time, type
install.packages('ProfileParsimony')

## Once the package has been installed, load it using
library(ProfileParsimony)

## Data can be read from a nexus file (note the restrictions detailed above):
my.data <- read.nexus.data('C:/path/to/filename.nex')

## Alternatively you can use a built-in dataset:
data(SigSut); my.data <- SigSut.data

## A contrast matrix translates the tokens used in your dataset to the character states to 
##    which they correspond: for example decoding 'A' to {01}.
##    For more details, see the 'phangorn-specials' vignette in the phangorn package, accesible 
##    by typing '?phangorn' in the R prompt and navigating to index > package vignettes.

contrast.matrix <- matrix(data=c(
# 0 1 -  # Each column corresponds to a character-state
  1,0,0, # Each row corresponds to a token, here 0, denoting the character-state set {0} 
  0,1,0, # 1 | {1}
  0,0,1, # - | {-}
  1,1,0, # A | {01}
  1,1,0, # + | {01}
  1,1,1  # ? | {01-}
), ncol=3, byrow=TRUE); # ncol should correspond to the number of columns in the matrix
dimnames(contrast.matrix) <- list(
  c(0, 1, '-', 'A', '+', '?'), # A list of the tokens corresponding to each rows in the contrast matrix
  c(0, 1, '-') # A list of the character-states corresponding to the columns in the contrast matrix
)

## To see the annotated contrast matrix, type
contrast.matrix

## Apply the contrast matrix, using the phyDat format...
my.phyDat <- PhyDat(my.data, type='USER', contrast=contrast.matrix)

## ... and optimize the data for analysis 
my.prepdata <- PrepareDataFitch(my.phyDat)

## Specify the names of the outgroup taxa
my.outgroup <- c('taxon1', 'taxon2')

## Load a bifurcating tree,
tree <- read.nexus('treename.nex')
## or generate a random starting tree, 
tree <- RandomTree(my.phyDat)
## or use neighbour joining to generate a starting tree
tree <- root(nj(dist.hamming(my.phyDat)), my.outgroup, resolve.root=TRUE); tree$edge.length <- NULL;

## View the starting tree by typing
plot(tree)

## Calculate the tree's parsimony score
FitchInfoFast(tree, my.prepdata)

## Search for a better tree
better.tree <- TreeSearch(tree, my.prepdata, my.outgroup, method='SPR')
## Try the parsimony ratchet (Nixon, 1999)
better.tree <- Pratchet(better.tree, my.prepdata, my.outgroup, maxhits=50, k=20)
## The default parameters may not be enough to find the most parsimonious tree; type 
##    ?Pratchet to view all search parameters.

## View the results
plot(better.tree)

## Collect all trees that are suboptimal by up to 1.5 bits
## The trick here is to use large values of k and pratchiter, and small values of searchhits
## and searchiter, so that many runs don't quite hit the optimal tree.
suboptimals <- Pratchet(RandomTree(my.prepdata), PrepareDataFitch(my.prepdata), 
                        all=TRUE, outgroup=my.outgroup, suboptimal=1.5, k=250, 
                        searchhits=10, searchiter=50, pratchiter=5000, rearrangements='TBR')

## Calculate and display a consensus tree that includes slightly suboptimal trees
plot(my.consensus <- consensus(suboptimals))
}
}
