\name{unreliability_MI}
\alias{unreliability_MI}
\title{Unreliability estimation and normalized mean intensities (MI) scores calculation 
  for Infinium probes on the Illumina MethylationEPIC microarray v1.0}

\description{
  Estimating probes 'unreliability' scores by simulating the beta-value
  distribution influenced by technical noise (which is getting from Negative Control Probes from Green and Red channels separately) on different levels of M (Methylated) and U (Unmethylated)
  intensities signals. High 'unreliability' values are associated with
  low reproducibility and low repeatability of beta-values and are due to low normalized mean intensities (MI) level
  (which is also calculated within a given function). Based on MI and unreliability scores, function also builds a suggestion for removing low-reliable probes.
}
\usage{
  unreliability_MI(RGset, samples, grid_max_intenisty = 5000, grid_step = 100, number_beta_generated = 1000)
}
\arguments{
  \item{RGset}{minfi RGChannelSet object}
  \item{samples}{(optional) subset of samples used to estimate all parameters (all samples by default)}
  \item{grid_max_intenisty}{(optional) maximum value of artificial intensity (M and U) of Reliability Map grid (5000 by default)} 
  \item{grid_step}{(optional) reliability Map grid step (of range from 0 to \code{grid_max_intenisty}) (100 by default)} 
  \item{number_beta_generated}{(optional) the number of generated beta-values (with noise) that are used to estimate the MAD of
    their distribution (1000 by default) for each M and U pair} 
}
\value{
  Return data frame of probes (\emph{probe} column) with indiacating their type (\emph{type_of_probe} column 
and \emph{channel} column for Type I probes); with unreliability score (\emph{unreliability} column), 
with MI (\emph{MI} column) and with recommendation for exculsion 
(\emph{recommended_for_exclusion} with values \emph{1} - if should be ecluded, \emph{0} - if not).
}
\references{
  Nazarenko T, Vavourakis CD, Jones A, Evans A, Watson AW, Brandt K,
  Carter C, Zaikin A, Herzog C, Widschwendter M.
  \emph{Technical and biological sources of unreliability of Infinium
  type II probes of the Illumina MethylationEPIC BeadChip microarray} (2023)
}
\author{
  Tatiana Nazarenko, Charlotte Vavourakis
}
\examples{
  \dontrun{
  # To run on all samples with default settings:
  out <- unreliability_MI(RGset)
}
  \dontrun{
  # To run on specific subset of \emph{samples} (samples_subset) with default settings:
  out <- unreliability_MI(RGset, samples = samples_subset)}
}
\keyword{unreliability}
\keyword{MI}
