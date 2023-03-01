load("./data/probesII_EPICv1.RData")
load("./data/probesI_EPICv1.RData")


unreliability_map_estimation <- function(noise_matrix, grid_max_intenisty, grid_step, number_beta_generated) {
  noise_green <- subset(noise_matrix, noise_matrix$type =="green")
  noise_red <- subset(noise_matrix, noise_matrix$type =="red")

  red_grid <- seq(0,grid_max_intenisty, grid_step)
  green_grid <- seq(0,grid_max_intenisty, grid_step)

  print("Unreliability Map calculation:")
  print("Beta generation...")

  beta_generator <- data.frame()

  for(i in 1:length(red_grid)) {
    for(j in 1:length(green_grid)) {
      red = red_grid[i]
      green = green_grid[j]
      noise_number <- sample(1:nrow(noise_green), number_beta_generated)
      noised_green = noise_green[noise_number,]$noise + green
      noised_red = noise_red[noise_number,]$noise + red
      beta_temp <- (noised_green)/(noised_green + noised_red)
      beta_generator <- rbind(beta_generator, data.frame(i,j,red, green,mean(noised_green),mean(noised_red), t(beta_temp)))
    }
  }
  print("Unreliability values calculation...")
  beta_generator <- data.frame(beta_generator)
  beta_generator_matrix <- beta_generator[7:(ncol(beta_generator))]

  range_beta_distribution <- apply(beta_generator_matrix,1,function(x) quantile(x, probs=c(0.025,0.975),na.rm=TRUE))
  range_beta_distribution <- data.frame(t(range_beta_distribution))

  unreliability_map <- cbind(beta_generator[,1:6], (range_beta_distribution))
  unreliability_map$q <- unreliability_map$X97.5. - unreliability_map$X2.5.
  unreliability_map$mean.noise_green. <- beta_generator$mean.noised_green.
  unreliability_map$mean.noise_red. <- beta_generator$mean.noised_red.
  unreliability_map$green <- beta_generator$green
  unreliability_map$red <- beta_generator$red
  unreliability_map$number = seq(1, nrow(unreliability_map))
  return(unreliability_map)
}



unreliability_calculation <- function(noise_matrix, samples,green_array, red_array,probes,unreliability_map, grid_max_intenisty, grid_step) {

  if (!requireNamespace("stringr", quietly = TRUE)) {
    stop(
      "Package stringr must be installed to use this function.",
      call. = FALSE
    )
  }

  mean_green_noise <- mean(subset(noise_matrix, noise_matrix$type =="green")$noise)
  mean_red_noise <- mean(subset(noise_matrix, noise_matrix$type =="red")$noise)

  probes_names <- as.vector(probes$probe)

  unreliability_array <-data.frame(cg = probes_names)

  n_steps = length(seq(0,grid_max_intenisty, grid_step))

  for(i in 1:length(samples)) {
    print(stringr::str_c("sample ",i,"/",length(samples)))

    sample <- samples[i]

    sample_probes_array <- probes
    sample_probes_array$green = as.vector(t(green_array[probes_names, sample]))
    sample_probes_array$red = as.vector(t(red_array[probes_names, sample]))

    sub_sample_probes_array <- subset(sample_probes_array, sample_probes_array$green < grid_max_intenisty & sample_probes_array$red < grid_max_intenisty)

    sample_red_adjusted_on_noise <- ifelse(sub_sample_probes_array$red - mean_red_noise <0,0,sub_sample_probes_array$red - mean_red_noise)
    sample_green_adjusted_on_noise <- ifelse(sub_sample_probes_array$green - mean_green_noise <0,0,sub_sample_probes_array$green - mean_green_noise)

    ureliability_map_closest_cell_number = (sample_red_adjusted_on_noise%/%grid_step)*n_steps + (sample_green_adjusted_on_noise%/%grid_step) + 1
    sample_probes_array$q <- 0
    sample_probes_array[rownames(s), ]$q <- unreliability_map[ureliability_map_closest_cell_number,]$q
    unreliability_array <- cbind(unreliability_array, sample_probes_array$q)

  }
  colnames(real_pat) <- c("cg",samples)
  rownames(unreliability_array) <- probes_names

  probes_unreliability<- apply(unreliability_array[2:ncol(unreliability_array)], 1, function(x) mean(x, na.rm=T))
  return(probes_unreliability)
}

unreliability_MI <- function(probes, RGset, noise_set="p", samples, grid_max_intenisty = 5000, grid_step = 100, number_beta_generated = 1000) {
  rownames(probes) <- probes$probe

  if (!requireNamespace("minfi", quietly = TRUE)) {
    stop(
      "Package minfi must be installed to use this function.",
      call. = FALSE
    )
  }

  MSet <- minfi::preprocessRaw(RGset)
  green_array <- minfi::getMeth(MSet)
  red_array <- minfi::getUnmeth(MSet)

  common_probes <- intersect(rownames(green_array), as.vector(probes$probe))
  probes <- probes[common_probes,]

  if(missing(samples)) {
  samples <- colnames(green_array)
}

  now <- Sys.time()
  print("Start...")
  if(noise_set == "Y") {
    noise_probes <- rownames(subset(probes, probes$CHR == "Y"))
  } else {
    if(noise_set == "p") {
      print("p-value calculating...")
      detP_rgset <- minfi::detectionP(RGset, type = "m+u")
      pvalue_array <- data.frame(detP_rgset)
      if(length(intersect(samples, colnames(pvalue_array))) != 0) {
        countp_001 <- apply(pvalue_array[,samples],1,function(x) sum(x>0.01)/length(x))
      }
      if(length(intersect(paste0("X",samples), colnames(pvalue_array))) != 0) {
        countp_001 <- apply(pvalue_array[,paste0("X",samples)],1,function(x) sum(x>0.01)/length(x))
      }
      noise_probes <- intersect(names(countp_001[countp_001>=0.5]), as.vector(probes$probe))
    } else {
          stop(
      "Please, select CORRECT method of noise probes selection: p -- by p-value; Y -- by Y chromosomes probes (applicable only for female samples)",
      call. = FALSE
    )
    }
  }


  if (length(intersect(samples, colnames(green_array))) == 0) {
    stop(
      "Loaded samples do not match by names with GREEN and RED arrays column names",
      call. = FALSE
    )}

  if (length(noise_probes) == 0) {
              stop(
      "0 noise probes were detect. Please, check you probes data frame (it should contain the maximum possible set of probes of the used array).",
      call. = FALSE
    )
    } else {
    print(str_c(length(noise_probes)," probes will be used for noise estimation"))
    }


  noise_probes <- intersect(noise_probes, as.vector(probes$probe))
  noise_probes_df <- probes[noise_probes,]
  noise_probes_df$mean_green <- apply(green_array[noise_probes, samples],1,function(x) mean(x, na.rm=T))
  noise_probes_df$mean_red <- apply(red_array[noise_probes, samples],1,function(x) mean(x, na.rm=T))

  noise_probes_common_intensity <- noise_probes_df$mean_green+noise_probes_df$mean_red
  names(noise_probes_common_intensity) <- noise_probes
  selected_noise_probes <- names(noise_probes_common_intensity[noise_probes_common_intensity<2*density(noise_probes_common_intensity)$x[which.max(density(noise_probes_common_intensity)$y)]] )

  green_noise <- green_array[selected_noise_probes, samples]
  red_noise <- red_array[selected_noise_probes, samples]

  noise_matrix <- data.frame()
  for(i in 1: length(selected_noise_probes)) {
    temp <- rbind(data.frame(type="green", noise = as.vector(t(green_noise[selected_noise_probes[i],]))),
                  data.frame(type="red", noise = as.vector(t(red_noise[selected_noise_probes[i],]))))
    noise_matrix <- rbind(noise_matrix, temp)
  }

  unreliability_map <- unreliability_map_estimation(noise_matrix, grid_max_intenisty, grid_step, number_beta_generated)
  print("Unreliability calculation for all data...")
  unreliability <- unreliability_calculation(noise_matrix, samples,green_array, red_array,probes,unreliability_map,grid_max_intenisty, grid_step)

  probes$unreliability <- unreliability
  print("MI calculating...")
  probes <- get_MI(probes, green_array, red_array, samples)
  print("Finish!")
  print(Sys.time() - now)
  return(probes)
}

get_MI <- function(probes, green_array, red_array, samples) {
  samples_df <- data.frame(samples)
  samples_df$mean_sample_green <- apply(green_array[as.vector(probes$probe),samples], 2, function(x) mean(x, na.rm=T))
  samples_df$mean_sample_red <- apply(red_array[as.vector(probes$probe),samples], 2, function(x) mean(x, na.rm=T))

  samples_df$mean_sample_intensity <- samples_df$mean_sample_green + samples_df$mean_sample_red

  probes_norm_green <- apply(green_array[as.vector(probes$probe),samples], 1, function(x) mean(x/samples_df$mean_sample_intensity))
  probes_norm_red <- apply(red_array[as.vector(probes$probe),samples], 1, function(x) mean(x/samples_df$mean_sample_intensity))

  probes$MI <- (probes_norm_green + probes_norm_red)

  return(probes)
}
