# --------------------------------------------------------------------------------------------------------
# Install data.table (most recent version)
# --------------------------------------------------------------------------------------------------------
# devtools::install_github('Rdatatable/data.table')
# --------------------------------------------------------------------------------------------------------
# Install h2o (most recent version)
# --------------------------------------------------------------------------------------------------------
# if ("package:h2o" %in% search()) detach("package:h2o", unload=TRUE)
# if ("h2o" %in% rownames(installed.packages())) remove.packages("h2o")
# # Next, download H2O package dependencies:
# pkgs <- c("methods","statmod","stats","graphics","RCurl","jsonlite","tools","utils")
# new.pkgs <- setdiff(pkgs, rownames(installed.packages()))
# if (length(new.pkgs)) install.packages(new.pkgs)
# # Download and install the H2O package for R:
# install.packages("h2o", type="source", repos=(c("https://s3.amazonaws.com/h2o-release/h2o/master/3636/R")))
# --------------------------------------------------------------------------------------------------------
# Install stremr
# --------------------------------------------------------------------------------------------------------
# devtools::install_github('osofr/stremr', build_vignettes = FALSE)

