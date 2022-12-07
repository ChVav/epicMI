
<!-- README.md is generated from README.Rmd. Please edit that file -->

# epicMI

<!-- badges: start -->
<!-- badges: end -->

Estimate unreliability and MI scores for Infinium type II probes on the
Illumina MethylationEPIC microarray v1.0, either using detP values or
probes on the Y chromosome.

## Requirements

R \>= 3.5.0 <br> minfi \>= 1.42.0 <br> stringr

## Installation

``` r
if(!require(devtools)) install.packages("devtools")
devtools::install_github("ChVav/epicMI")
```

## Example

``` r
library(epicMI)
out <- main_unreliability_MI(probesII_EPICv1, RGset)
```
