\name{unreliability_MI}
\alias{unreliability_MI}
\title{Unreliability estimation and normalized mean intensities (MI) scores calculation 
  for Infinium type II probes on the Illumina MethylationEPIC microarray v1.0}

\description{
  Method for estimating probes 'unreliability' scores by simulating the beta-value
  distribution influenced by technical noise on different levels of green and red
  intensities signals. High "unreliability" values are associated with
  low reproducibility and low repeatability of beta-values and are due to low MI level
  (which is also calculated within a given function).
}
\usage{
  unreliability_MI(probesII, RGset, noise_set,samples)
}
\arguments{
  \item{probes}{A dataframe of type II or type I probes (of specific version of Illumina
    Methylation BeadChip array) which should have at least 2 column:
    - 'probe' column (name of probes);
    - 'CHR' column (with chromosome annotation).
    
    In package embedded probesII_EPICv1 dataframe (for type II probes of Illumina
    MethylationEPIC BeadChip v1.0) and probesI_EPICv1 dataframe (for type I probes
    of Illumina MethylationEPIC BeadChip v1.0). These dataframes also have additional
    column 'removed_to_EPICv2', where "1" - indicated that probe of Illumina
    MethylationEPIC BeadChip v1.0 was removed to Illumina MethylationEPIC BeadChip v2.0,
    and "0" - was not.}
  
  \item{RGset}{minfi RGChannelSet object}
  
  \item{noise_set}{Method to select noise probes,
    either 'p' (default) - select noise probes by p-value (if >= 50 percent of samples have p-value > 0.01)
    or 'Y' (only for samples from female subjects) - select noise probes as Y chromosomes.}
  
  \item{samples}{Optional, subset of samples used to estimate noise. Useful for using
    Y-noise method (for pre-selection of female samples). If missing, then all RGset IDs are
    used.}
   
  \item{grid_max_intenisty}{Maximum value of intensity (green and red) of Reliability Map grid (5000 by default).} 
    
  \item{grid_step}{Reliability Map grid step (of range from 0 to \code{grid_max_intenisty}) (100 by default).} 
    
  \item{number_beta_generated}{The number of generated beta-values (with noise) that are used to estimate the range of
    their distribution (1000 by default).} 
  
}
\value{
  Original probes dataframe with 2 columns ('unreliability' and 'MI') added.
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
  ## run calculation on probes type II of Illumina MethylationEPIC BeadChip v1.0 (probesII_EPICv1)
  ## with noise probes detected by p_values (default) on ALL samples (default) 
  ## with grid_max_intenisty = 5000 (default), grid_step = 100 (default), number_beta_generated = 1000 (default)
  \dontrun{out <- unreliability_MI(probesII_EPICv1, RGset)}
  
  ## run calculation with noise probes as Ychr probes on female samples
  \dontrun{out <- unreliability_MI(probesII_EPICv1, RGset, noise_set="Y", samples = female_samples)}
  
  ## run calculation with noise probes as Ychr probes on female samples
  ## without probesII type which were removed for EPICv2
  \dontrun{out <- unreliability_MI(subset(probesII_EPICv1,probesII_EPICv1$removed_to_EPICv2 == 0) , RGset, noise_set="Y", samples = female_samples)}
}
\keyword{unreliability}
\keyword{MI}
