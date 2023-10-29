
<!-- README.md is generated from README.Rmd. Please edit that file -->

# epicMI

<!-- badges: start -->
<!-- badges: end -->

Using the *unreliability_MI* function, calculate normalized mean intensities (MI), estimate unreliability and scores and identify threshold for removing Infinium probes on the
Illumina MethylationEPIC microarray v1.0.

## Requirements

R \>= 3.5.0 <br> minfi \>= 1.42.0 <br> stringr <br> dplyr <br> ggpubr <br> ggplot2


## Installation

``` r
if(!require(devtools)) install.packages("devtools")
devtools::install_github("ChVav/epicMI")
```

## Example
Loading package and preparing RGset:

``` r
library(epicMI)
RGset <- read.metharray.exp(targets)
```

``` r
out <- unreliability_MI(RGset, samples, grid_max_intenisty = 5000, grid_step = 100, number_beta_generated = 1000)
```

By default, *unreliability_MI* function calculate normalized mean intensities (MI), estimate unreliability and scores and identify threshold for removing probes estimate unreliability and calculate MI scores on all samples, using *p-noise* method, with estimation on Reliability Map grid: M(0,*'grid_max_intenisty'*) x U(0,*'grid_max_intenisty'*) (where by default *'grid_max_intenisty'* = 5000), 
with *'grid_step'* = 100 and *'number_beta_generated'* = 1000:

``` r
out <- unreliability_MI(RGset, samples)
```
*Note* Method can be apply to other Illumina Methylation microarrays (450k, EPICv2.0 etc)



