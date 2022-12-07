\name{probesII_EPICv1}
\alias{probesII_EPICv1}
\docType{data}
\title{
Basic annotation of Infinium Type II probes on Illumina MethylationEPICv1.0 microarray.
}
\description{
Basic annotation of Infinium Type II probes on Illumina MethylationEPICv1.0 microarray. For these probes the MI and unreliability scores will be calculated by default by unreliability_MI_typeII().
}
\usage{data("probesII_EPICv1")}
\format{
  A data frame with 723722 observations on the following 9 variables.
  \describe{
    \item{\code{probe}}{probe name}
    \item{\code{countC}}{number of C in probe sequence}
    \item{\code{countA}}{number of A in probe sequence}
    \item{\code{countR}}{number of purines in probe sequence}
    \item{\code{CHR}}{targetted chromosome, sourced from original manifest}
    \item{\code{chrXY}}{"X" or "Y" sex chromosome; "other" autosome, sourced from original manifest}
    \item{\code{snps}}{SNP in position "0", "1", "both", or not a SNP in position 0 or 1 from target CpG, sourced from original manifest}
    \item{\code{Region}}{CpG context, sourced from original manifest}
    \item{\code{removed_to_EPICv2}}{"1" probe not on EPIC v2.0, "0" probe also on EPIC v2.0}
  }
}
\references{
https://support.illumina.com/downloads/infinium-methylationepic-v1-0-product-files.html (Mar 13, 2020)
}
\keyword{datasets}
