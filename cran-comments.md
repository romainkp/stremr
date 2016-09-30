## Initial submission to CRAN:

## Test environments:
* local OS X install, R 3.2.4
* Ubuntu 12.04 R 3.3.1
* (offline -- could not be tested) win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* Package suggested but not available for checking: ‘h2oEnsemble’

This package provides an important added functionality for stremr and is available for easy installation from github. A reference to its github repository has been added to the field "Additional_repositories:". Appropriate error checks and messages have been implemented throughout the package. 
Currently, the error message instructs the user on how to install it directly from github, i.e.,
"Package h2oEnsemble is needed for modeling with SuperLearner.
Please install it by typing this into R terminal:
  library(devtools)
  install_github(\"h2oai/h2o-3/h2o-r/ensemble/h2oEnsemble-package\")"
Finally, the maintainer of h2oEnsemble has promised to release the package to CRAN very soon.

