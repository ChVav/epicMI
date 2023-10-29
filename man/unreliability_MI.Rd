\name{unreliability_MI}
\alias{unreliability_MI}
\title{Unreliability estimation and normalized mean intensities (MI) scores calculation 
  for Infinium probes on the Illumina MethylationEPIC microarray v1.0}

\description{
  Estimating probes 'unreliability' scores by simulating the beta-value
  distribution influenced by technical noise (which is getting from Negative Control Probes from Green and Red channels separately) on different levels of M (Methylated) and U (Unmethylated)
  intensities signals. High "unreliability" values are associated with
  low reproducibility and low repeatability of beta-values and are due to low normalized mean intensities (MI) level
  (which is also calculated within a given function). Based on MI and unreliability scores, function also builds a suggestion for removing low-reliable probes.

}
\usage{
  unreliability_MI(RGset, samples, grid_max_intenisty = 5000, grid_step = 100, number_beta_generated = 1000)
}
\arguments{
  \item{RGset}{minfi RGChannelSet object}
  \item{samples}{(optional) subset of samples used to estimate all parameters (all samples by default)}
  
  * parameters related to the construction of the Reliability Map *
  \item{grid_max_intenisty}{(optional) maximum value of artificial intensity (M and U) of Reliability Map grid (5000 by default).} 
    
  \item{grid_step}{(optional) reliability Map grid step (of range from 0 to \code{grid_max_intenisty}) (100 by default).} 
    
  \item{number_beta_generated}{(optional) the number of generated beta-values (with noise) that are used to estimate the MAD of
    their distribution (1000 by default) for each M and U pair.} 
  
  \item{RGset}{minfi RGChannelSet object} 
}
\value{
  Return data frame of probes (\emph{probe} column) with indiacating their type (\emph{type_of_probe} column and \emph{channel} column for Type I probes);
with unreliability score (\emph{unreliability} column), with MI (\emph{MI} column) and with recommendation for exculsion 
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
  ## To run on all samples with default settings:

  \dontrun{out <- unreliability_MI(RGset)}

  ## To run on specific subset of \code{samples} ('samples_subset') with default settings:

  \dontrun{out <- unreliability_MI(RGset, samples = samples_subset)}

  ## To run with specific values for \code{grid_max_intenisty} or \code{grid_step] or \code{number_beta_generated} ('specific_grid_max_intenisty', 'specific_grid_step' or 'specific_number_beta_generated')

  \dontrun{out <- unreliability_MI(RGset, grid_max_intenisty = specific_grid_max_intenisty, grid_step = specific_grid_step, number_beta_generated = specific_number_beta_generated)}
  
}
\keyword{unreliability}
\keyword{MI}
