
# load remotes package to install older package versions
library(remotes)

# set an existing directory in which the packages can be installed
package_dir <- "~/my_directory/"

# specify empty subdirectories where to install the package versions
registrOld_dir <- paste0(package_dir, "registr_1.0.0/")
registrNew_dir <- paste0(package_dir, "registr_2.1.4/")
fdasrvf_dir    <- paste0(package_dir, "fdasrvf_1.9.4/")



# package installations ---------------------------------------------------
### 1) registr 1.0.0, based on lme4 1.1-23, gamm4 0.2-6
# registr 1.0.0 is saved locally since the original 1.0.0 version was
# edited to include the current version of the convergence criterion
# to be comparable to registr 2.0.7.
# Mainly, this change comprises the use of a convergence threshold of 0.0001
# instead of 0.01 in the main while loop of registr::register_fpca().
remotes::install_local("registr_1.0.0_newConvergenceCriterion/", lib = registrOld_dir)
remotes::install_version("lme4",  version = "1.1-23", lib = registrOld_dir)
remotes::install_version("gamm4", version = "0.2-6",  lib = registrOld_dir)


### 2) registr 2.1.4, based on lme4 1.1.26, gamm4 0.2.7
remotes::install_github("julia-wrobel/registr@gamm4", lib = registrNew_dir)
remotes::install_version("lme4", version = "1.1.26",  lib = registrNew_dir)
remotes::install_github("r-gam/gamm4@issue-1",        lib = registrNew_dir)


### 3) fdasrvf x.x.x
remotes::install_version("fdasrvf", version = "1.9.4", lib = fdasrvf_dir)



# Test to load / unload the package versions ------------------------------
### load old registr version
library("lme4",    lib.loc = registrOld_dir)
library("gamm4",   lib.loc = registrOld_dir)
library("registr", lib.loc = registrOld_dir)
packageVersion("lme4")
packageVersion("gamm4")
packageVersion("registr")

# unload the old registr version
detach("package:registr", unload = TRUE)
detach("package:gamm4",   unload = TRUE)
detach("package:lme4",    unload = TRUE)


### load new registr version
library("lme4",    lib.loc = registrNew_dir)
library("gamm4",   lib.loc = registrNew_dir)
library("registr", lib.loc = registrNew_dir)
packageVersion("lme4")
packageVersion("gamm4")
packageVersion("registr")


### load fdasrvf
library("fdasrvf", lib.loc = fdasrvf_dir)
packageVersion("fdasrvf")
