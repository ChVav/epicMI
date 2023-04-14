load("./data/probesII_EPICv1.RData")
load("./data/probesI_EPICv1.RData")

#' @title
#' Unreliability estimation and normalized mean intensities (MI) scores calculation
#' for Infinium probes on the Illumina MethylationEPIC microarray v1.0
#'
#' @aliases unreliability_MI
#'
#' @description
#' Method for estimating probes 'unreliability' scores by simulating the beta-value
#' distribution influenced by technical noise on different levels of green and red
#' intensities signals. High "unreliability" values are associated with
#' low reproducibility and low repeatability of beta-values and are due to low MI level
#' (which is also calculated within a given function).
#'
#' @param probes
#' A dataframe of type II or type I probes (of specific version of Illumina
#' Methylation BeadChip array) which should have at least 2 column:
#' - 'probe' column (name of probes);
#' - 'CHR' column (with chromosome annotation).
#'
#' In package embedded probesII_EPICv1 dataframe (for type II probes of Illumina
#' MethylationEPIC BeadChip v1.0) and probesI_EPICv1 dataframe (for type I probes
#' of Illumina MethylationEPIC BeadChip v1.0). These dataframes also have additional
#' column 'removed_to_EPICv2', where "1" - indicated that probe of Illumina
#' MethylationEPIC BeadChip v1.0 was removed to Illumina MethylationEPIC BeadChip v2.0,
#' and "0" - was not.
#'
#' @param RGset
#' An object of class \code{RGChannelSet}.
#'
#' @param noise_set
#' Choice of a method of noise probes selection ('p' = default, 'Y'). 'p'-noise method -
#' probes, which are failed (p-value > 0.01) on 50% samples, where p-values calculated
#' by \code{detectionP} function (\code{minfi} package); Y-noise method - probes on the
#' Y chromosome (for using only on female samples).
#'
#' @param samples
#' Optional parameter of sample list (should be identical to IDs in RGset). Useful for using
#' Y-noise method (for pre-selection of female samples). If missing, then all RGset IDs are
#' used as a sample list.
#'
#' @param grid_max_intenisty
#' Maximum value of intensity (green and red) of Reliability Map grid (5000 by default).
#'
#' @param grid_step
#' Reliability Map grid step (of range from 0 to \code{grid_max_intenisty}) (100 by default).
#'
#' @param number_beta_generated
#' The number of generated beta-values (with noise) that are used to estimate the range of
#' their distribution (1000 by default).
#'
#' @return \code{probes} dataframe with 2 columns ('unreliability' and 'MI') added.
#'
#' @references
#' Nazarenko T, Vavourakis CD, Jones A, Evans A, Watson AW, Brandt K,
#' Carter C, Zaikin A, Herzog C, Widschwendter M.
#' \emph{Technical and biological sources of unreliability of Infinium
#' type II probes of the Illumina MethylationEPIC BeadChip microarray} (2023)
#'
#'
#' @examples
#' out <- unreliability_MI(probesII_EPICv1, RGset)
#'
#' out <- unreliability_MI(probesII_EPICv1, RGset, noise_set="Y", samples)
#'
#' @export unreliability_MI
#'
unreliability_MI <- function(RGset, samples, list_of_noise_probes, grid_max_intenisty = 5000, grid_step = 100, number_beta_generated = 1000) {


  if (!requireNamespace("minfi", quietly = TRUE)) stop("Package minfi must be installed to use this function.", call. = FALSE)

  now <- Sys.time()
  print("Start...")

  probesII <- minfi::getProbeInfo(RGset, type = 'II')$Name

  probesI_green <- minfi::getProbeInfo(RGset, type = 'I-Green')$Name
  probesI_red<- minfi::getProbeInfo(RGset, type = 'I-Red')$Name

  probes <- rbind(data.frame(probe = probesI_green, channel = "green", type_of_probe = "I"),
                  data.frame(probe = probesI_red, channel = "red", type_of_probe = "I"),
                  data.frame(probe = probesII, channel = "", type_of_probe = "II"))
  rownames(probes) <- probes$probe

  MSet <- minfi::preprocessRaw(RGset)
  green_array <- minfi::getMeth(MSet)
  red_array <- minfi::getUnmeth(MSet)

  if(missing(samples)) {
    samples <- colnames(green_array)
  } else {
    if (length(intersect(samples, colnames(green_array))) == 0) stop("Loaded samples do not match by names with GREEN and RED arrays column names", call. = FALSE)
  }

  if(missing(list_of_noise_probes)) {
    print("p-value calculating...")
    detP_rgset <- minfi::detectionP(RGset, type = "m+u")
    pvalue_array <- data.frame(detP_rgset)
    if(length(intersect(samples, colnames(pvalue_array))) != 0) {
      countp_001 <- apply(pvalue_array[, samples], 1, function(x) sum(x > 0.01) / length(x))
    }
    if(length(intersect(paste0("X", samples), colnames(pvalue_array))) != 0) {
      countp_001 <- apply(pvalue_array[, paste0("X", samples)], 1, function(x) sum(x > 0.01) / length(x))
    }
    noise_probes <- names(countp_001[countp_001 >= 0.5])
  } else {
    noise_probes <- list_of_noise_probes
    if (length(noise_probes) == 0) {
      stop("0 noise probes were detected. Please, check your list_of_noise_probes", call. = FALSE)
    }
  }
  noise_probes_df <- probes[noise_probes,]

  green_noise <- green_array[c(subset(noise_probes_df, noise_probes_df$type_of_probe == "II")$probe,
                               subset(noise_probes_df, noise_probes_df$type_of_probe == "I" & noise_probes_df$channel == "green")$probe), samples]
  red_noise <- red_array[c(subset(noise_probes_df, noise_probes_df$type_of_probe == "II")$probe,
                               subset(noise_probes_df, noise_probes_df$type_of_probe == "I" & noise_probes_df$channel == "red")$probe), samples]

  green_noise_df <- data.frame(x=as.numeric(unlist(green_noise)))
  red_noise_df <- data.frame(x=as.numeric(unlist(red_noise)))

  print(green_noise_df$x)
  selected_green_noise <- green_noise_df$x[green_noise_df$x < 2 * density(green_noise_df$x)$x[which.max(density(green_noise_df$x)$y)]]
  selected_red_noise <- red_noise_df$x[red_noise_df$x < 2 * density(red_noise_df$x)$x[which.max(density(red_noise_df$x)$y)]]

  noise_matrix <- rbind(data.frame(type = "green", noise = selected_green_noise),
                        data.frame(type = "red", noise = selected_red_noise))

  results_on_typeI_green <- specific_type_probe_calculation(noise_matrix, "I-green", samples, grid_max_intenisty, grid_step, number_beta_generated, probesI_green, green_array, red_array)
  results_on_typeI_red <- specific_type_probe_calculation(noise_matrix, "I-red", samples, grid_max_intenisty, grid_step, number_beta_generated, probesI_red, green_array, red_array)
  results_on_typeII <- specific_type_probe_calculation(noise_matrix, "II", samples, grid_max_intenisty, grid_step, number_beta_generated, probesII, green_array, red_array)

  probes <- rbind(results_on_typeI_green, results_on_typeI_red,results_on_typeII)
  print("*Finish*")
  print(str_c("Time: ", Sys.time() - now))
  return(probes)
}

specific_type_probe_calculation <- function(noise_matrix, type_of_probes, samples, grid_max_intenisty, grid_step, number_beta_generated, probe_names, green_array, red_array) {
  if(type_of_probes == "I-green") {
    probes <- data.frame(probe = probe_names, channel = "green", type_of_probe = "I")
    rownames(probes) <- probes$probe
  }
  if(type_of_probes == "I-red") {
    probes <- data.frame(probe = probe_names, channel = "red", type_of_probe = "I")
    rownames(probes) <- probes$probe
  }
  if(type_of_probes == "II") {
    probes <- data.frame(probe = probesII, channel = "", type_of_probe = "II")
    rownames(probes) <- probesII
  }

  print(str_c("Calculation for type ", type_of_probes,"..."))

  unreliability_map <- unreliability_map_estimation(noise_matrix, type_of_probes, grid_max_intenisty, grid_step, number_beta_generated)
  print(str_c("...type ",type_of_probes," probes unreliability calculation ..."))
  unreliability <- unreliability_calculation(noise_matrix, type_of_probes, samples, green_array, red_array, probes, unreliability_map, grid_max_intenisty, grid_step)
  probes$unreliability <- unreliability

  print(str_c("...type ",type_of_probes, " probes MI calculating..."))
  probes <- get_MI(probes, green_array, red_array, samples)
  return(probes)
}

unreliability_map_estimation <- function(noise_matrix, type_of_probes, grid_max_intenisty, grid_step, number_beta_generated) {

  if(type_of_probes == "I-green") {
    noise_green <- subset(noise_matrix, noise_matrix$type == "green")
    noise_red <- subset(noise_matrix, noise_matrix$type == "green")
  }
  if(type_of_probes == "I-red") {
    noise_green <- subset(noise_matrix, noise_matrix$type == "red")
    noise_red <- subset(noise_matrix, noise_matrix$type == "red")
  }
  if(type_of_probes == "II") {
    noise_green <- subset(noise_matrix, noise_matrix$type == "green")
    noise_red <- subset(noise_matrix, noise_matrix$type == "red")
  }

  red_grid <- seq(0, grid_max_intenisty, grid_step)
  green_grid <- seq(0, grid_max_intenisty, grid_step)

  print("Unreliability Map calculation:")
  print("Beta generation...")
  beta_generator <- data.frame()
  for(i in 1:length(red_grid)) {
    for(j in 1:length(green_grid)) {
      red = red_grid[i]
      green = green_grid[j]
      noise_number_green <- sample(1:nrow(noise_green), number_beta_generated)
      noise_number_red <- sample(1:nrow(noise_red), number_beta_generated)
      noised_green = noise_green[noise_number_green, ]$noise + green
      noised_red = noise_red[noise_number_red, ]$noise + red
      beta_temp <- (noised_green) / (noised_green + noised_red)
      beta_generator <- rbind(beta_generator, data.frame(i, j, red, green, mean(noised_green), mean(noised_red), t(beta_temp)))
    }
  }

  print("Unreliability values calculation...")
  beta_generator <- data.frame(beta_generator)
  beta_generator_matrix <- beta_generator[7:(ncol(beta_generator))]
  range_beta_distribution <- apply(beta_generator_matrix, 1, function(x) quantile(x, probs=c(0.025, 0.975), na.rm = TRUE))
  range_beta_distribution <- data.frame(t(range_beta_distribution))

  unreliability_map <- cbind(beta_generator[, 1:6], (range_beta_distribution))
  unreliability_map$q <- unreliability_map$X97.5. - unreliability_map$X2.5.
  unreliability_map$mean.noise_green. <- beta_generator$mean.noised_green.
  unreliability_map$mean.noise_red. <- beta_generator$mean.noised_red.
  unreliability_map$green <- beta_generator$green
  unreliability_map$red <- beta_generator$red
  unreliability_map$number = seq(1, nrow(unreliability_map))
  return(unreliability_map)
}

unreliability_calculation <- function(noise_matrix, type_of_probes, samples, green_array, red_array, probes, unreliability_map, grid_max_intenisty, grid_step) {

  if (!requireNamespace("stringr", quietly = TRUE)) stop("Package stringr must be installed to use this function.", call. = FALSE)

  if(type_of_probes == "I-green") {
    mean_green_noise <- mean(subset(noise_matrix, noise_matrix$type == "green")$noise)
    mean_red_noise <- mean(subset(noise_matrix, noise_matrix$type == "green")$noise)
  }
  if(type_of_probes == "I-red") {
    mean_green_noise <- mean(subset(noise_matrix, noise_matrix$type == "red")$noise)
    mean_red_noise <- mean(subset(noise_matrix, noise_matrix$type == "red")$noise)
  }
  if(type_of_probes == "II") {
    mean_green_noise <- mean(subset(noise_matrix, noise_matrix$type == "green")$noise)
    mean_red_noise <- mean(subset(noise_matrix, noise_matrix$type == "red")$noise)
  }

  probes_names <- as.vector(probes$probe)
  unreliability_array <-data.frame(cg = probes_names)

  n_steps = length(seq(0, grid_max_intenisty, grid_step))

  for(i in 1:length(samples)) {
    print(stringr::str_c("sample ", i, "/", length(samples)))

    sample <- samples[i]

    sample_probes_array <- probes
    sample_probes_array$green = as.vector(t(green_array[probes_names, sample]))
    sample_probes_array$red = as.vector(t(red_array[probes_names, sample]))

    sub_sample_probes_array <- subset(sample_probes_array, sample_probes_array$green < grid_max_intenisty & sample_probes_array$red < grid_max_intenisty)

    sample_red_adjusted_on_noise <- ifelse(sub_sample_probes_array$red - mean_red_noise < 0, 0, sub_sample_probes_array$red - mean_red_noise)
    sample_green_adjusted_on_noise <- ifelse(sub_sample_probes_array$green - mean_green_noise < 0, 0, sub_sample_probes_array$green - mean_green_noise)

    reliability_map_closest_cell_number = (sample_red_adjusted_on_noise %/% grid_step)*n_steps + (sample_green_adjusted_on_noise %/% grid_step) + 1
    sample_probes_array$q <- 0
    sample_probes_array[rownames(sub_sample_probes_array), ]$q <- unreliability_map[reliability_map_closest_cell_number,]$q
    unreliability_array <- cbind(unreliability_array, sample_probes_array$q)

  }
  colnames(unreliability_array) <- c("cg", samples)
  rownames(unreliability_array) <- probes_names

  probes_unreliability <- apply(unreliability_array[2:ncol(unreliability_array)], 1, function(x) mean(x, na.rm=T))
  return(probes_unreliability)
}

get_MI <- function(probes, green_array, red_array, samples) {
  samples_df <- data.frame(samples)
  samples_df$mean_sample_green <- apply(green_array[as.vector(probes$probe), samples], 2, function(x) mean(x, na.rm = T))
  samples_df$mean_sample_red <- apply(red_array[as.vector(probes$probe), samples], 2, function(x) mean(x, na.rm = T))

  samples_df$mean_sample_intensity <- samples_df$mean_sample_green + samples_df$mean_sample_red

  probes_norm_green <- apply(green_array[as.vector(probes$probe), samples], 1, function(x) mean(x / samples_df$mean_sample_intensity))
  probes_norm_red <- apply(red_array[as.vector(probes$probe), samples], 1, function(x) mean(x / samples_df$mean_sample_intensity))

  probes$MI <- (probes_norm_green + probes_norm_red)

  return(probes)
}
