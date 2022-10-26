# KGENP -------------------------------------------------------------------

RNP <- function(sim, obs){

  # Calculate mean sim and obs
  mean.sim = mean(sim, na.rm=TRUE)
  mean.obs = mean(obs, na.rm=TRUE)
  # Calculate normalized flow duration curves
  fdc.sim = sort(sim / (mean.sim * length(sim)))
  fdc.obs = sort(obs / (mean.obs * length(obs)))
  # Calculate alpha component
  RNP.alpha = 1 - 0.5 * sum(abs(fdc.sim - fdc.obs))
  # Calculate beta component
  RNP.beta = mean.sim / mean.obs
  # Calculate r component
  RNP.r = cor(sim, obs, method="spearman")
  # Calculate KGENP
  KGENP <-  1 - sqrt((RNP.alpha - 1)^2 + (RNP.beta - 1)^2 + (RNP.r - 1)^2)
  # Create list output
  KGENP_output <- list("KGENP" = KGENP,
                       "r" = RNP.r,
                       "Beta" = RNP.beta,
                       "Alpha" = RNP.alpha)
  
  # Return Non-Parametric Efficiency value
  return(KGENP_output)
}

# KGE'' -------------------------------------------------------------------

KGE_m2 <- function(sim, obs, out.type = "single") {
  
  # Calculate mean sim and obs
  mean.sim <- mean(sim, na.rm=TRUE)
  mean.obs <- mean(obs, na.rm=TRUE)
  # Calculate sd sim and obs
  sd.sim <- sd(sim, na.rm=TRUE)
  sd.obs <- sd(obs, na.rm=TRUE)
  # Calculate r component
  r <- cor(sim, obs, method = "pearson")
  # Calculate alpha component
  alpha <-  sd.sim / sd.obs
  # Calculate beta component
  beta <- (mean.sim - mean.obs) / sd.obs
  # Calculate KGE''
  KGE_m2_score <- 1 - sqrt(beta ^ 2 + (alpha - 1) ^ 2 + (r - 1) ^ 2)
  
  # Output
  if (out.type == "full") {
    KGE_m2_output <- list("KGE_m2" = KGE_m2_score,
                          "r" = r,
                          "Beta_m" = beta,
                          "Alpha" = alpha)
    
    return(KGE_m2_output)
  }
  
  return(KGE_m2_score)
}

# LME ---------------------------------------------------------------------

LME <- function(sim, obs, out.type = "single") {
  
  # Calculate mean sim and obs
  mean.sim <- mean(sim, na.rm=TRUE)
  mean.obs <- mean(obs, na.rm=TRUE)
  # Calculate sd sim and obs
  sd.sim <- sd(sim, na.rm=TRUE)
  sd.obs <- sd(obs, na.rm=TRUE)
  # Calculate r component
  r <- cor(sim, obs, method = "pearson")
  # Calculate alpha component
  alpha <-  sd.sim / sd.obs
  # Calculate beta component
  beta <- mean.sim / mean.obs
  # Calculate LME
  LME_score <- 1 - sqrt(((r * alpha - 1) ^ 2 + (beta - 1) ^ 2))
  
  # Output
  if (out.type == "full") {
    LME_output <- list("LME" = LME_score,
                       "Alpha" = alpha,
                       "r" = r,
                       "r_alpha" = r * alpha,
                       "Beta" = beta)
    
    return(LME_output)
  }
  
  return(LME_score)
}


# LCE ---------------------------------------------------------------------

LCE <- function(sim, obs, out.type = "single") {
  
  # Calculate mean sim and obs
  mean.sim <- mean(sim, na.rm=TRUE)
  mean.obs <- mean(obs, na.rm=TRUE)
  # Calculate sd sim and obs
  sd.sim <- sd(sim, na.rm=TRUE)
  sd.obs <- sd(obs, na.rm=TRUE)
  # Calculate r component
  r <- cor(sim, obs, method = "pearson")
  # Calculate alpha component
  alpha <-  sd.sim / sd.obs
  # Calculate beta component
  beta <- mean.sim / mean.obs
  # Calculate LCE
  LCE_score <- 1 - sqrt(((r * alpha - 1) ^ 2 + (r / alpha - 1) ^ 2 + (beta - 1) ^ 2))
  
  # Output
  if (out.type == "full") {
    LCE_output <- list("LCE" = LCE_score,
                       "Alpha" = alpha,
                       "r" = r,
                       "r_alpha" = r * alpha,
                       "r_on_alpha" = r / alpha,
                       "Beta" = beta)
    
    return(LCE_output)
  }
  
  return(LCE_score)
}


# NSE ---------------------------------------------------------------------

NSE_custom <- function(sim, obs) {
  
  # Calculate mean sim and obs
  mean.sim <- mean(sim, na.rm=TRUE)
  mean.obs <- mean(obs, na.rm=TRUE)
  # Calculate sd sim and obs
  sd.sim <- sd(sim, na.rm=TRUE)
  sd.obs <- sd(obs, na.rm=TRUE)
  # Calculate alpha component
  alpha <-  sd.sim / sd.obs
  # Calculate beta component
  beta <- (mean.sim - mean.obs) / sd.obs
  # Calculate r component
  r <- cor(sim, obs, method = "pearson")
  # Calculate NSE
  NSE_score <- 2 * alpha * r - alpha ^ 2 - beta ^ 2
  
  # Output
  NSE_output <- list("NSE" = NSE_score,
                     "r" = r,
                     "2_alpha_r" = 2 * alpha * r,
                     "Beta_m" = beta,
                     "Alpha" = alpha)
  
  return(NSE_output)
}

# KGE_sf ---------------------------------------------------------------------

KGE_sf <- function(sim, obs, sf = c("cor" = 1, "var" = 1, "bias" = 1), method = "2009") {
  
  # Get scaling factors
  cor_sf <- sf[["cor"]]
  var_sf <- sf[["var"]]
  bias_sf <- sf[["bias"]]
  # Mean values
  mean.sim <- mean(sim, na.rm=TRUE)
  mean.obs <- mean(obs, na.rm=TRUE)
  # Standard deviations
  sigma.sim <- sd(sim, na.rm=TRUE)
  sigma.obs <- sd(obs, na.rm=TRUE)
  # Pearson product-moment correlation coefficient
  r <- cor(sim, obs, method = "pearson") 
  # Alpha is a measure of relative variability between simulated and observed values (See Ref1)
  Alpha <- sigma.sim / sigma.obs
  # Calculate beta component
  beta <- mean.sim / mean.obs
  # CV.sim is the coefficient of variation of the simulated values [dimensionless]
  # CV.obs is the coefficient of variation of the observations [dimensionless]
  CV.sim <- sigma.sim / mean.sim
  CV.obs <- sigma.obs / mean.obs
  # Gamma is the variability ratio, which is used instead of Alpha
  Gamma <- CV.sim / CV.obs
  
  # Calculate KGE_sf
  if (method == "2012") {
    KGE_sf <- 1 - sqrt((bias_sf * (beta - 1))^2 + 
                         (var_sf * (Gamma - 1))^2 + 
                         (cor_sf * (r - 1))^2) 
  } else if (method == "2009") {
    KGE_sf <- 1 - sqrt((bias_sf * (beta - 1))^2 + 
                         (var_sf * (Alpha - 1))^2 + 
                         (cor_sf * (r - 1))^2) 
  }
  
  # Output
  return(KGE_sf)
}