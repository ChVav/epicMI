load("./data/probesII_EPICv1.RData")
load("./data/probesI_EPICv1.RData")


beta_generator <- function(Noise_matrix) {
  noiseG <- subset(Noise_matrix, Noise_matrix$type =="green")
  noiseR <- subset(Noise_matrix, Noise_matrix$type =="red")
  
  N = 1000
  red_grid <- seq(0,5000, 100)
  green_grid <- seq(0,5000, 100)
  beta_generator <- data.frame()
  green_generator <- data.frame()
  red_generator <- data.frame()

  print("Unreliability Map calculation:")
  print("Beta generation...")
  for(i in 1:length(red_grid)) {
    for(j in 1:length(green_grid)) {
      red = red_grid[i]
      green = green_grid[j]
      noise_number <- sample(1:nrow(noiseG), N)
      noised_green = noiseG[noise_number,]$noise + green
      noised_red = noiseR[noise_number,]$noise + red
      beta_temp <- (noised_green)/(noised_green + noised_red)
      beta_generator <- rbind(beta_generator, data.frame(i,j,red, green,mean(noised_green),mean(noised_red), t(beta_temp)))
    }
  }
  print("Unreliability values calculation...")
  beta_generator <- data.frame(beta_generator)
  colnames(beta_generator)
  beta_generator2 <- beta_generator[7:(ncol(beta_generator))]
  
  qqq <- apply(beta_generator2,1,function(x) quantile(x, probs=c(0.025,0.975),na.rm=TRUE))
  qqq <- data.frame(t(qqq))
  beta_generator2 <- cbind(beta_generator[,1:6], (qqq))
  beta_generator2$q <- beta_generator2$X97.5. - beta_generator2$X2.5.
  beta_generator2$mean.noise_green. <- beta_generator$mean.noised_green.
  beta_generator2$mean.noise_red. <- beta_generator$mean.noised_red.
  beta_generator2$green <- beta_generator$green
  beta_generator2$red <- beta_generator$red
  beta_generator2$number = seq(1, nrow(beta_generator2))
  return(beta_generator2)
}



unreliability_calculation <- function(Noise_matrix, samples,GREEN, RED,probes,beta_generator2) {
  
  if (!requireNamespace("stringr", quietly = TRUE)) {
    stop(
      "Package stringr must be installed to use this function.",
      call. = FALSE
    )
  }
  

  
  noiseG <- subset(Noise_matrix, Noise_matrix$type =="green")
  noiseR <- subset(Noise_matrix, Noise_matrix$type =="red")
  
  tII <- as.vector(probes$probe)
  meanR_noise <- mean(noiseR$noise)
  meanG_noise <- mean(noiseG$noise)
  
  real_pat<-data.frame(cg = tII)
  
  Nsteps = 51
  step=100
  
  temp_probes <- probes
  
  for(i in 1:length(samples)) {
    print(stringr::str_c("sample ",i,"/",length(samples)))
    start_time <- Sys.time()
    
    ssample <- samples[i]
    
    temp_probes$green = as.vector(t(GREEN[tII, ssample]))
    temp_probes$red = as.vector(t(RED[tII, ssample]))
    
    s <- subset(temp_probes, temp_probes$green < 5000 & temp_probes$red < 5000)
    green <- s$green
    red <- s$red
    
    red <- ifelse(red - meanR_noise <0,0,red - meanR_noise)
    green <- ifelse(green - meanG_noise <0,0,green - meanG_noise)
    k = (red%/%step)*Nsteps + (green%/%step) + 1
    temp_probes$k <- 0
    temp_probes[rownames(s), ]$k <- beta_generator2[k,]$q
    real_pat <- cbind(real_pat, temp_probes$k)
    end_time <- Sys.time()
    #print(end_time - start_time)
  }
  colnames(real_pat) <- c("cg",samples)
  rownames(real_pat) <- tII
  
  mrel<- apply(real_pat[2:ncol(real_pat)], 1, function(x) mean(x, na.rm=T))
  return(mrel[tII])
}

unreliability_MI <- function(probes, RGset, noise_set,samples) {
  
  if (!requireNamespace("minfi", quietly = TRUE)) {
    stop(
      "Package minfi must be installed to use this function.",
      call. = FALSE
    )
  }
  
  MSet <- minfi::preprocessRaw(RGset)
  GREEN <- minfi::getMeth(MSet)
  RED <- minfi::getUnmeth(MSet)
  
  common_probes <- intersect(rownames(GREEN), as.vector(probes$probe))
  probes <- probes[common_probes,]

  if(missing(samples)) {
  samples <- colnames(GREEN)
} 
  
  if(missing(noise_set)) {
    noise_set <- "p"
  }
  
  now <- Sys.time()
  print("Start...")
  if(noise_set == "Y") {
    NOISE_PROBES <- rownames(subset(probes, probes$CHR == "Y"))
  } else {
    if(noise_set == "p") {
      print("p-value calculating...")
      detP_rgset <- minfi::detectionP(RGset, type = "m+u")
      PVAL_V2 <- data.frame(detP_rgset)
      countp_001 <- apply(PVAL_V2[,paste0("X",samples)],1,function(x) sum(x>0.01)/length(x))
      NOISE_PROBES <- intersect(names(countp_001[countp_001>=0.5]), as.vector(probes$probe))
    } else {
          stop(
      "ERROR: Please, select CORRECT method of noise probes selection: p -- by p-value; Y -- by Y chromosomes probes (applicable only for female samples)",
      call. = FALSE
    )
    }
  }
  print("!!!")                        
  print(length(intersect(samples, colnames(GREEN)))
                          
  if (length(intersect(samples, colnames(GREEN))) == 0) {
    stop(
      str_c("ERROR: loaded samples (",samples[0],", ",samples[1],", ",samples[2],", ..." ,") do not match by names 
with GREEN and RED arrays columns (",colnames(GREEN)[0],",v",colnames(GREEN)[1],",v",colnames(GREEN)[2]," ..." ,")"),
      call. = FALSE
    )} 
        
  if (length(NOISE_PROBES) == 0) {
              stop(
      "ERROR: 0 noise probes were detect. Please, check you probes data frame (it should contain the maximum possible set of probes of the used array.)",
      call. = FALSE
    )
    } else {
    print(str_c(length(NOISE_PROBES)," probes will be used for noise estimation")) 
    }
                   

  NOISE_PROBES <- intersect(NOISE_PROBES, as.vector(probes$probe))
  NOISE_PROBES_df <- probes[NOISE_PROBES,]
  NOISE_PROBES_df$mean_green <- apply(GREEN[NOISE_PROBES, samples],1,function(x) mean(x, na.rm=T))
  NOISE_PROBES_df$mean_red <- apply(RED[NOISE_PROBES, samples],1,function(x) mean(x, na.rm=T))

  v_NOISE_PROBES <- NOISE_PROBES_df$mean_green+NOISE_PROBES_df$mean_red
  names(v_NOISE_PROBES) <- NOISE_PROBES
  good_NOISE_PROBES <- names(v_NOISE_PROBES[v_NOISE_PROBES<2*density(v_NOISE_PROBES)$x[which.max(density(v_NOISE_PROBES)$y)]] )
  
  GREEN_NOISE <- GREEN[good_NOISE_PROBES, samples]
  RED_NOISE <- RED[good_NOISE_PROBES, samples]
  
  Noise_matrix <- data.frame()
  for(i in 1: length(good_NOISE_PROBES)) {
    temp <- rbind(data.frame(type="green", noise = as.vector(t(GREEN_NOISE[good_NOISE_PROBES[i],]))),
                  data.frame(type="red", noise = as.vector(t(RED_NOISE[good_NOISE_PROBES[i],]))))
    Noise_matrix <- rbind(Noise_matrix, temp)
  }
  
  bg <- beta_generator(Noise_matrix)
  print("unreliability calculation for all data...")
  uc <- unreliability_calculation(Noise_matrix, samples,GREEN, RED,probes,bg)
  
  probes$unreliability <- uc
  print("MI calculating...")
  probes <- get_MI(probes, GREEN, RED, samples)
  print("FINISH!")
  print(Sys.time() - now)
  return(probes)
}

get_MI <- function(probes, GREEN, RED, samples) {
  samples_df <- data.frame(samples)
  samples_df$meanG <- apply(GREEN[as.vector(probes$probe),samples], 2, function(x) mean(x, na.rm=T))
  samples_df$meanR <- apply(RED[as.vector(probes$probe),samples], 2, function(x) mean(x, na.rm=T))
  
  samples_df$meanI <- samples_df$meanG + samples_df$meanR

  green_mean_v1_norm <- apply(GREEN[as.vector(probes$probe),samples], 1, function(x) mean(x/samples_df$meanI))
  red_mean_v1_norm <- apply(RED[as.vector(probes$probe),samples], 1, function(x) mean(x/samples_df$meanI))

  probes$MI <- (green_mean_v1_norm + red_mean_v1_norm)
  
  return(probes)
}
