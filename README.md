
<!-- README.md is generated from README.Rmd. Please edit that file -->

# epicMI

<!-- badges: start -->
<!-- badges: end -->

The function evaluates unreliability and calculates MI values based on "noise" probes, which can be selected in two ways:

Estimate unreliability and calculate MI scores for Infinium type II probes on the
Illumina MethylationEPIC microarray v1.0, based on "noise" probes, which can be selected in two ways:
- *p-noise* method - probes, which are failed (p-value > 0.01) on 50% samples, where p-values calculated by *detectionP* function (**minfi** package) 
- *Y-noise* method - probes on the Y chromosome (for using only on female samples).

## Requirements

R \>= 3.5.0 <br> minfi \>= 1.42.0 <br> stringr

## Installation

``` r
if(!require(devtools)) install.packages("devtools")
devtools::install_github("ChVav/epicMI")
```


## Example

### *p-noise* method 
By default, *unreliability_MI* function estimate unreliability and calculate MI scores on all samples, using *p-noise* method:
``` r
library(epicMI)
out <- unreliability_MI(probesII_EPICv1, RGset)
```
In this example we are using probes II Type of EPIC microarray v1.0: data frame **probesII_EPICv1**, which will loaded with package.

For using own set of probes (for example on for estimation of different version Illumina Methylation array) this set should have the same structure and at least have *'probe'* column (name of probes) and *'CHR'* column (with chromosome annotation).
Additionaly, we added column *'removed_to_EPICv2'* where "1" - indicated that probe of EPIC microarray v1.0 was removed to EPIC microarray v2.0, and "0" - was not.

Alternatively, this method can be run as:
``` r
out <- unreliability_MI(probesII_EPICv1, RGset, noise_set="p")
```
or on selected samples as
``` r
out <- unreliability_MI(probesII_EPICv1, RGset, noise_set="p", samples),
```
where samples names should be the same as in RGset.

### *Y-noise* method 

For using *Y-noise* method here are 2 options to run *unreliability_MI* function:
``` r
out <- unreliability_MI(probesII_EPICv1, RGset, noise_set="Y"),
```
if all the samples are female, 
or 
``` r
out <- unreliability_MI(probesII_EPICv1, RGset, noise_set="Y", samples),
```
where samples names should be the same as in RGset and should be only female samples.
