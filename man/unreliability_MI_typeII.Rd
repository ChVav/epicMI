\name{unreliability_MI_typeII}
\alias{unreliability_MI_typeII}
\title{unreliability_MI_typeII
}
\description{
Calculates the mean intensity (MI) and unreliability scores for Infininium type II probes on the Illumina MethylationEPIC array v1.0.
}
\usage{
unreliability_MI_typeII(probesII, RGset, noise_set,samples)
}
\arguments{
  \item{probesII}{Default = probesII_EPICv1; or user-defined set of probes.
  Should have 'probe' column (with name of CpG) and CHR column (with chromosome annotation)}
  \item{RGset}{minfi RGChannelSet object}
  \item{noise_set}{Method to select Type II noise probes,
  either 'p' (default) - select noise probes by p-value (if >= 50 percent of samples have p-value > 0.01)
  or 'Y' (only for samples from female subjects) - select noise probes as Y chromosomes.}
  \item{samples}{Optional, subset of samples used to estimate noise.}
}
\value{
Outputs a dataframe with for each Infinium type II probe the estimated MI and unreliability scores.
}
\references{
Please cite xxx.
}
\author{
Tatiana Nazarenko, Charlotte Vavourakis
}
\examples{
## run calculation with noise probes detected by p_values (default) on ALL samples (default)\
\dontrun{out <- main_unreliability_MI(probesII_EPICv1, RGset)}
## or
\dontrun{out <- main_unreliability_MI(probesII_EPICv1, RGset, noise_set="p")}

## run calculation with noise probes as Ychr probes on female samples
\dontrun{out <- main_unreliability_MI(probesII_EPICv1, RGset, noise_set="Y", samples)}

## run calculation with noise probes as Ychr probes on female samples
## without probesII type which were removed for EPICv2
\dontrun{out <- main_unreliability_MI(subset(probesII_EPICv1,probesII_EPICv1$removed_to_EPICv2 == 0) , RGset, noise_set="Y", samples)}
}
\keyword{unreliability}
\keyword{MI}
