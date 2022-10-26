#' Calculate performance criteria between sim and obs values.
#'
#' @param sim Numerical vector. A series of simulated values.
#' @param obs Numerical vector. A series of observed values.
#'
#' @return A list with 9 performance criteria and their components
#' KGE, KGE', KGE'', KGEnp, DE, LME, LCE, NSE, dr
#' @export
#'
#' @examples

cperf <- function(sim, obs) {
  
  # Prepare output
  perfc <- list()
  
  ## KGE
  KGE <- hydroGOF::KGE(sim, obs, out.type = "full")
  
  perfc$KGE$score <- KGE$KGE.value
  perfc$KGE$beta <- KGE$KGE.elements[["Beta"]]
  perfc$KGE$alpha <- KGE$KGE.elements[["Alpha"]]
  perfc$KGE$r <- KGE$KGE.elements[["r"]]
  
  ## KGE modified (KGE')
  KGE_m <- hydroGOF::KGE(sim, obs, out.type = "full", method = "2012")
  
  perfc$KGE_m$score <- KGE_m$KGE.value
  perfc$KGE_m$beta <- KGE_m$KGE.elements[["Beta"]]
  perfc$KGE_m$gamma <- KGE_m$KGE.elements[["Gamma"]]
  perfc$KGE_m$r <- KGE_m$KGE.elements[["r"]]
  
  ## KGE_modified2 (KGE'')
  KGE_m2 <- KGE_m2(sim, obs, out.type = "full")
  
  perfc$KGE_m2$score <- KGE_m2$KGE_m2
  perfc$KGE_m2$beta_n <- KGE_m2$Beta_m
  perfc$KGE_m2$alpha <- KGE_m2$Alpha
  perfc$KGE_m2$r <- KGE_m2$r
  
  ## non-parametric KGE (KGENP)
  KGENP <- RNP(sim, obs)
  
  perfc$KGENP$score <- KGENP$KGENP
  perfc$KGENP$beta <- KGENP$Beta
  perfc$KGENP$alpha_np <- KGENP$Alpha
  perfc$KGENP$r <- KGENP$r
  
  ## LME
  LME <- LME(sim, obs, out.type = "full")
  
  perfc$LME$score <- LME$LME
  perfc$LME$beta <- LME$Beta
  perfc$LME$r_alpha <- LME$r_alpha
  
  ## LCE
  LCE <- LCE(sim, obs, out.type = "full")
  
  perfc$LCE$score <- LCE$LCE
  perfc$LCE$beta <- LCE$Beta
  perfc$LCE$r_alpha <- LCE$r_alpha
  perfc$LCE$r_on_alpha <- LCE$r_on_alpha
  
  ## NSE
  NSE <- NSE_custom(sim, obs)
  
  perfc$NSE$score <- NSE$NSE
  perfc$NSE$beta_n <- NSE$Beta_m
  perfc$NSE$alpha <- NSE$Alpha
  perfc$NSE$`2_alpha_r` <- NSE$`2_alpha_r`
  
  ## dr
  dr <- hydroErr$dr(sim, obs)
  
  perfc$dr$score <- dr
  
  ## DE
  brel <- de$calc_brel_mean(obs, sim, sort=sort)
  barea <- de$calc_bias_area(de$calc_brel_res(obs, sim, sort=sort))
  r <- de$calc_temp_cor(obs, sim)
  de <- 1 - de$calc_de(obs, sim)
  
  DE <- list("DE" = de,
             "r" = r,
             "brel" = brel,
             "barea" = barea)
  
  perfc$DE$score <- DE$DE
  perfc$DE$brel <- DE$brel
  perfc$DE$barea <- DE$barea
  perfc$DE$r <- DE$r
  
  # Output
  return(perfc)
}
