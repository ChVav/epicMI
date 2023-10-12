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

unreliability_MI <- function(RGset, samples, grid_max_intenisty = 5000, grid_step = 100, number_beta_generated = 1000) {
  
  # Check if the minfi package is available
  if (!requireNamespace("minfi", quietly = TRUE)) stop("Package minfi must be installed to use this function.", call. = FALSE)
  
  # Record the starting time for execution
  now <- Sys.time()
  print("Start...")
  
  # Extract probe names for different types
  probesII <- minfi::getProbeInfo(RGset, type = 'II')$Name
  probesI_green <- minfi::getProbeInfo(RGset, type = 'I-Green')$Name
  probesI_red<- minfi::getProbeInfo(RGset, type = 'I-Red')$Name
  
  # Combine probe information into a structured dataframe
  probes <- rbind(data.frame(probe = probesI_green, channel = "green", type_of_probe = "I"),
                  data.frame(probe = probesI_red, channel = "red", type_of_probe = "I"),
                  data.frame(probe = probesII, channel = "", type_of_probe = "II"))
  rownames(probes) <- probes$probe
  
  # Preprocess the raw intensity data
  MSet <- minfi::preprocessRaw(RGset)
  green_array <- minfi::getMeth(MSet)
  red_array <- minfi::getUnmeth(MSet)
  
  # Check if samples are provided, if not, use all available columns
  if(missing(samples)) {
    samples <- colnames(green_array)
  } else {
    if (length(intersect(samples, colnames(green_array))) == 0) stop("Loaded samples do not match by names with GREEN and RED arrays column names", call. = FALSE)
  }
  
  print("Noise from negative control probes extraction...")
  
  # Extract noise from negative control probes
  probes_control <- minfi::getProbeInfo(RGset, type = 'Control')
  probes_control <- data.frame(probes_control[2:nrow(probes_control),])
  probes_control_negative <- subset(probes_control,probes_control$Type == "NEGATIVE")
  
  # Extract values from Green and Red channels
  GREEN_channel <- getGreen(RGset)
  RED_channel<- getRed(RGset)
  
  # Intersection of noise probes in channels
  green_n <- intersect(rownames(GREEN_channel), probes_control_negative$Address)
  red_n <-  intersect(rownames(RED_channel), probes_control_negative$Address)
  
  # Create data frames for green and red noise
  green_noise_df_nc <- data.frame(x=as.numeric(unlist(GREEN_channel[green_n,])))
  red_noise_df_nc <- data.frame(x=as.numeric(unlist(RED_channel[red_n,])))
  
  # Select noise values below 99th percentile
  selected_green_noise <- green_noise_df_nc$x[green_noise_df_nc$x < quantile(green_noise_df_nc$x, probs = c(0.99))]
  selected_red_noise <- red_noise_df_nc$x[red_noise_df_nc$x < quantile(red_noise_df_nc$x, probs = c(0.99))]
  
  # Create a noise matrix
  noise_matrix <- rbind(data.frame(type = "green", noise = selected_green_noise),
                        data.frame(type = "red", noise = selected_red_noise))
  
  # Calculate unreliability and MI for different probe types
  results_on_typeI_green <- specific_type_probe_calculation(noise_matrix, "I-green", samples, grid_max_intenisty, grid_step, number_beta_generated, probesI_green, green_array, red_array)
  results_on_typeI_red <- specific_type_probe_calculation(noise_matrix, "I-red", samples, grid_max_intenisty, grid_step, number_beta_generated, probesI_red, green_array, red_array)
  results_on_typeII <- specific_type_probe_calculation(noise_matrix, "II", samples, grid_max_intenisty, grid_step, number_beta_generated, probesII, green_array, red_array)
  
  # Calculate speed change point for different probe types
  speed_change_point_out_Igreen <- speed_change_point(results_on_typeI_green, "I_green")
  results_on_typeI_green <- speed_change_point_out_Igreen[[1]]
  fig_typeI_green <- speed_change_point_out_Igreen[[2]]
  
  speed_change_point_out_Ired <- speed_change_point(results_on_typeI_red, "I_red")
  results_on_typeI_red <- speed_change_point_out_Ired[[1]]
  fig_typeI_red <- speed_change_point_out_Ired[[2]]
  
  speed_change_point_out_II <- speed_change_point(results_on_typeII, "II")
  results_on_typeII <- speed_change_point_out_II[[1]]
  fig_typeII <- speed_change_point_out_II[[2]]
  
  # Combine the results for different probe types
  probes <- rbind(results_on_typeI_green, results_on_typeI_red,results_on_typeII)
  figure <- ggarrange(fig_typeI_green, fig_typeI_red,fig_typeII, nrow=3)
  plot(figure)
  print("*Finish*")
  print(str_c("Time: ", Sys.time() - now))
  return(probes)
}
??ggarrange
# Function to calculate unreliability and mutual information for a specific type of probes
specific_type_probe_calculation <- function(noise_matrix, type_of_probes, samples, grid_max_intenisty, grid_step, number_beta_generated, probe_names, green_array, red_array) {
  
  # Create a data frame to store probe information based on the probe type
  if(type_of_probes == "I-green") {
    probes <- data.frame(probe = probe_names, channel = "green", type_of_probe = "I")
    rownames(probes) <- probes$probe
  }
  if(type_of_probes == "I-red") {
    probes <- data.frame(probe = probe_names, channel = "red", type_of_probe = "I")
    rownames(probes) <- probes$probe
  }
  if(type_of_probes == "II") {
    probes <- data.frame(probe = probe_names, channel = "", type_of_probe = "II")
    rownames(probes) <- probes$probe
  }
  
  print(str_c("Calculation for type ", type_of_probes,"..."))
  
  # Estimate the unreliability map for the specific probe type
  unreliability_map <- unreliability_map_estimation(noise_matrix, type_of_probes, grid_max_intenisty, grid_step, number_beta_generated)
  
  print(str_c("...type ",type_of_probes," probes unreliability calculation ..."))
  
  # Calculate unreliability for the specific probe type
  unreliability <- unreliability_calculation(noise_matrix, type_of_probes, samples, green_array, red_array, probes, unreliability_map, grid_max_intenisty, grid_step)
  probes$unreliability <- unreliability
  
  print(str_c("...type ",type_of_probes, " probes MI calculating..."))
  
  # Calculate MI for the specific probe type
  probes <- get_MI(probes, green_array, red_array, samples)
  return(probes)
}

unreliability_map_estimation <- function(noise_matrix, type_of_probes, grid_max_intenisty, grid_step, number_beta_generated) {
  
  # Check the type of probe and filter noise matrix accordingly
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
  
  # Create a grid for intensity values (red and green)
  red_grid <- seq(0, grid_max_intenisty, grid_step)
  green_grid <- seq(0, grid_max_intenisty, grid_step)
  
  print("Unreliability Map calculation:")
  print("Beta generation...")
  
  # Initialize a dataframe to store generated beta values
  beta_generator <- data.frame()
  for(i in 1:length(red_grid)) {
    for(j in 1:length(green_grid)) {
      red = red_grid[i]
      green = green_grid[j]
      
      # Get random noisy values
      noise_number_green <- sample(1:nrow(noise_green), number_beta_generated)
      noise_number_red <- sample(1:nrow(noise_red), number_beta_generated)
      
      # Generate noisy intensity values
      noised_green = noise_green[noise_number_green, ]$noise + green
      noised_red = noise_red[noise_number_red, ]$noise + red
      
      # Calculate generated beta values
      beta_temp <- (noised_green) / (noised_green + noised_red)
      
      # Append beta values to the beta generator dataframe
      beta_generator <- rbind(beta_generator, data.frame(i, j, red, green, mean(noised_green), mean(noised_red), t(beta_temp)))
    }
  }

  print("Unreliability values calculation...")
  beta_generator <- data.frame(beta_generator)
  beta_generator_matrix <- beta_generator[7:(ncol(beta_generator))]

  unreliability_map <- cbind(beta_generator[, 1:6])
  
  # Calculate unreliability values (median absolute deviation)
  unreliability_map$q <- unlist(apply(beta_generator_matrix, 1, function(x) mad(x, na.rm = TRUE)))
  
  # Add other necessary columns to the unreliability map
  unreliability_map$mean.noise_green. <- beta_generator$mean.noised_green.
  unreliability_map$mean.noise_red. <- beta_generator$mean.noised_red.
  unreliability_map$green <- beta_generator$green
  unreliability_map$red <- beta_generator$red
  unreliability_map$number = seq(1, nrow(unreliability_map))
  return(unreliability_map)
}

unreliability_calculation <- function(noise_matrix, type_of_probes, samples, green_array, red_array, probes, unreliability_map, grid_max_intenisty, grid_step) {
  
  # Check if the stringr package is available
  if (!requireNamespace("stringr", quietly = TRUE)) stop("Package stringr must be installed to use this function.", call. = FALSE)
  
  # Check the type of probe and filter noise matrix accordingly
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
  
  # Initialize a dataframe to store unreliability values
  unreliability_array <-data.frame(cg = probes_names)
  
  n_steps = length(seq(0, grid_max_intenisty, grid_step))
  
  # Iterate over samples to calculate unreliability for each probe
  for(i in 1:length(samples)) {
    sample = samples[i]
    print(stringr::str_c("sample ", i, "/", length(samples)))
    
    sample_probes_array <- probes
    
    sample_probes_array$green = as.vector(green_array[probes_names, sample])
    sample_probes_array$red = as.vector(red_array[probes_names, sample])
    
    sub_sample_probes_array <- subset(sample_probes_array, sample_probes_array$green < grid_max_intenisty & sample_probes_array$red < grid_max_intenisty)
    
    # Adjust sample intensity values based on noise
    sample_red_adjusted_on_noise <- ifelse(sub_sample_probes_array$red - mean_red_noise < 0, 0, sub_sample_probes_array$red - mean_red_noise)
    sample_green_adjusted_on_noise <- ifelse(sub_sample_probes_array$green - mean_green_noise < 0, 0, sub_sample_probes_array$green - mean_green_noise)
    
    # Determine the closest cell in the unreliability map
    reliability_map_closest_cell_number = (sample_red_adjusted_on_noise %/% grid_step)*n_steps + (sample_green_adjusted_on_noise %/% grid_step) + 1
    sample_probes_array$q <- 0
    sample_probes_array[rownames(sub_sample_probes_array), ]$q <- unreliability_map[reliability_map_closest_cell_number,]$q
    unreliability_array <- cbind(unreliability_array, sample_probes_array$q)
    
  }
  colnames(unreliability_array) <- c("cg", samples)
  rownames(unreliability_array) <- probes_names
  
  # Calculate the average unreliability for each probe
  probes_unreliability <- apply(unreliability_array[2:ncol(unreliability_array)], 1, function(x) mean(x, na.rm=T))
  return(probes_unreliability)
}

get_MI <- function(probes, green_array, red_array, samples) {
  
  # Create a dataframe to store sample information
  samples_df <- data.frame(samples)
  samples_df$mean_sample_green <- apply(green_array[as.vector(probes$probe), samples], 2, function(x) mean(x, na.rm = T))
  samples_df$mean_sample_red <- apply(red_array[as.vector(probes$probe), samples], 2, function(x) mean(x, na.rm = T))
  
  samples_df$mean_sample_intensity <- samples_df$mean_sample_green + samples_df$mean_sample_red
  
  # Normalize green and red intensities based on the mean sample intensity
  probes_norm_green <- apply(green_array[as.vector(probes$probe), samples], 1, function(x) mean(x / samples_df$mean_sample_intensity))
  probes_norm_red <- apply(red_array[as.vector(probes$probe), samples], 1, function(x) mean(x / samples_df$mean_sample_intensity))
  
  probes$MI <- (probes_norm_green + probes_norm_red)
  
  return(probes)
}

speed_change_point <- function(probes, type_of_probes) {
  # Build smoothed plot
  plot <- ggplot(probes, aes(MI, unreliability)) + geom_smooth(n = 200)  
  
  # Exctract smoothed curve
  points_df <- extract_smoothed_points(plot)
  
  # gradient calculation
  gradients <- diff(points_df$unreliability) / diff(points_df$MI)
  points_df$grad <- c(-1, gradients)
  
  # second gradient calculation
  gradients2 <- diff(points_df$grad) / diff(points_df$MI)
  points_df$grad2 <- c(-1, gradients2)
  
  # exclude first 2 strings (without exact gradients values)
  points_df <- points_df[(3:nrow(points_df)),]

  # find the point of maximum speed change
  points_df$grad_minus <- -points_df$grad2 
  df_max <- mutate(points_df, local.minima = if_else(lag(grad_minus) > grad_minus & lead(grad_minus) > grad_minus, TRUE, FALSE))
  first_max <- subset(df_max,df_max$local.minima == "TRUE")[1,]
  
  plot_final <- plot + geom_vline(xintercept = first_max$MI) + ggtitle(paste0("Thresh for ", type_of_probes, " type: MI = ", round(first_max$MI,4)))
  
  probes$recommended_for_exclusion <- ifelse(probes$MI > round(first_max$MI,4), 0, 1)
  return(list(probes, plot_final))
}
                           
extract_smoothed_points <- function(plot) {
  # Ensure that the plot has been created with ggplot2 and geom_smooth
  if (!inherits(plot, "ggplot")) {
    stop("Input must be a ggplot object with geom_smooth.")
  }
  
  # Build the plot to access its data
  plot_data <- ggplot_build(plot)$data
  
  # Extract the data for the smoothed curve
  smoothed_data <- plot_data[[1]]
  
  # Return the X and Y coordinates of the smoothed curve
  return(data.frame(MI = smoothed_data$x, unreliability = smoothed_data$y))
  
}


