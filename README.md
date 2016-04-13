# ProfileParsimony

ProfileParsimony is an R package, built on Phangorn, to implement the Profile Parsimony approach
to character weighting introduced by Faith and Trueman (2001).  It also provides an implementation
of Successive Approximations weighting (Farris 1969).

It incorporates modifications to phangorn that improve the speed of phylogenetic analysis, 
and adds support for TBR rearrangements.

The package will soon be compiled and uploaded to the CRAN repository.  
Meanwhile, you can install the latest version of the package into R thus:

```r
# Install the devtools package from CRAN
install.packages('devtools')

# Install the inapplicable package from github
devtools::install_github('ms609/ProfileParsimony')

# Load the package into R
library('ProfileParsimony')
```

## References
D. P. Faith, J. W. H. Trueman, Towards an inclusive philosophy for phylogenetic inference.
Syst. Biol. 50, 331–350 (2001).  <doi:10.1080/10635150118627>

Farris, J. S. (1969). A successive approximations approach to character weighting. 
Systematic Biology, 18(4), 374–385. <doi:10.2307/2412182>

K. P. Schliep, phangorn: phylogenetic analysis in R. Bioinformatics. 27, 592–593 (2011).
<doi:10.1093/bioinformatics/btq706>
