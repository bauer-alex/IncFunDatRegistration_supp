# Overview

This repository contains the code and data supplement for the paper
"Registration for Incomplete Non-Gaussian Functional Data"
(currently under revision).

The repository is structured as follows:

- `code` contains the code to reproduce all analyses and figures
- `data` contains the seismic datasets
- `figures` contains all figures for the main paper and the online appendix
- `IncRegHelpers` is an R package containing helper functions for the analyses
- `results` contains the results for the simulation study and the estimated model
object for the seismic application

After having cloned the repository to your local system, the `IncRegHelpers`
package can be installed using the following call:
```r
devtools::install("IncRegHelpers")
```