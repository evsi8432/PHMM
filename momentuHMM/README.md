# momentuHMM  [![R-CMD-check](https://github.com/bmcclintock/momentuHMM/workflows/R-CMD-check/badge.svg)](https://github.com/bmcclintock/momentuHMM/actions) [![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [![CRAN_Downloads](https://cranlogs.r-pkg.org/badges/momentuHMM)](https://cran.r-project.org/package=momentuHMM) [![CRAN_Downloads](https://cranlogs.r-pkg.org/badges/grand-total/momentuHMM)](https://cran.r-project.org/package=momentuHMM)

R package for maximum likelihood analysis of animal movement behavior using multivariate Hidden Markov Models 

Get started with the vignette: [Guide to using momentuHMM](https://cran.r-project.org/package=momentuHMM/vignettes/momentuHMM.pdf)

## Installation instructions

### CRAN release
The package is available at [![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/momentuHMM)](https://cran.r-project.org/package=momentuHMM). To install it:
``` R
install.packages("momentuHMM")
```

### Install from Github
To install the latest (stable) version of the package from Github: [![R-CMD-check](https://github.com/bmcclintock/momentuHMM/workflows/R-CMD-check/badge.svg)](https://github.com/bmcclintock/momentuHMM/actions)
``` R
library(remotes)
install_github("evsi8432/momentuHMM")
```

To install the latest (**unstable**) version of the package from Github: [![R-CMD-check](https://github.com/bmcclintock/momentuHMM/actions/workflows/r.yml/badge.svg?branch=develop)](https://github.com/bmcclintock/momentuHMM/actions/workflows/r.yml)
``` R
library(remotes)
install_github("evsi8432/momentuHMM@develop")
```

## References
McClintock, B.T., Michelot, T. (2018) [momentuHMM: R package for generalized hidden Markov models of animal movement](http://dx.doi.org/10.1111/2041-210X.12995). *Methods in Ecology and Evolution*, 9(6), 1518-1530.

McClintock, B.T., King R., Thomas L., Matthiopoulos J., McConnell B.J., Morales J.M. (2012) [A general discrete-time modeling framework for animal movement using multistate random walks](https://doi.org/10.1890/11-0326.1). *Ecological Monographs*, 82(3), 335-349.

McClintock, B.T. (2017) [Incorporating telemetry error into hidden Markov models of animal movement using multiple imputation](https://link.springer.com/article/10.1007/s13253-017-0285-6). *Journal of Agricultural, Biological, and Environmental Statistics*, 22(3), 249-269.
