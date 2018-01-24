[![Build Status](https://travis-ci.org/ms609/ProfileParsimony.svg?branch=master)](https://travis-ci.org/ms609/ProfileParsimony)
[![codecov](https://codecov.io/gh/ms609/ProfileParsimony/branch/master/graph/badge.svg)](https://codecov.io/gh/ms609/ProfileParsimony)
[![CRAN Status Badge](http://www.r-pkg.org/badges/version/ProfileParsimony)](https://cran.r-project.org/package=ProfileParsimony)
<!--[![CRAN Downloads](http://cranlogs.r-pkg.org/badges/ProfileParsimony)](https://cran.r-project.org/package=ProfileParsimony)-->
<!--[![Research software impact](http://depsy.org/api/package/cran/ProfileParsimony/badge.svg)](http://depsy.org/package/r/ProfileParsimony)-->
[![Project Status: Inactive – The project has reached a stable, usable state but is no longer being actively developed; support/maintenance will be provided as time allows.](http://www.repostatus.org/badges/latest/inactive.svg)](http://www.repostatus.org/#inactive)


# ProfileParsimony

ProfileParsimony was an R package, built on Phangorn, to implement the Profile Parsimony approach
to character weighting introduced by Faith and Trueman (2001).  It also provides an implementation
of Successive Approximations weighting (Farris 1969).

It incorporates modifications to phangorn that improve the speed of phylogenetic analysis, 
and adds support for TBR rearrangements.

It has now been integrated into the more capable package [TreeSearch](https://github.com/ms609/TreeSearch).


It is still possible to install the standalone package thus:
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
