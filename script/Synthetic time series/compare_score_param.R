compare_score_param <- function(obs, bad_bad, bad_good) {
  
  # Calculate score and get details of parameters for each transformation ----------------
  
  # KGE
  KGE_BB <- KGE(bad_bad, obs, out.type = "full")
  KGE_BG <- KGE(bad_good, obs, out.type = "full")
  
  # KGE'
  KGE_m_BB <- KGE(bad_bad, obs, out.type = "full", method = "2012")
  KGE_m_BG <- KGE(bad_good, obs, out.type = "full", method = "2012")
  
  # non-parametric KGE
  KGENP_BB <- RNP(bad_bad, obs)
  KGENP_BG <- RNP(bad_good, obs)
  
  # KGE''
  KGE_m2_BB <- KGE_m2(bad_bad, obs, out.type = "full")
  KGE_m2_BG <- KGE_m2(bad_good, obs, out.type = "full")
  
  # LME
  LME_BB <- LME(bad_bad, obs, out.type = "full")
  LME_BG <- LME(bad_good, obs, out.type = "full")
  
  # LCE
  LCE_BB <- LCE(bad_bad, obs, out.type = "full")
  LCE_BG <- LCE(bad_good, obs, out.type = "full")
  
  # NSE
  NSE_BB <- NSE_custom(bad_bad, obs)
  NSE_BG <- NSE_custom(bad_good, obs)
  
  # dr
  dr_BB <- hydroErr$dr(bad_bad, obs)
  dr_BG <- hydroErr$dr(bad_good, obs)
  
  # DE'
  brel_bb <- de$calc_brel_mean(obs, bad_bad, sort=sort)
  barea_bb <- de$calc_bias_area(de$calc_brel_res(obs, bad_bad, sort=sort))
  r_bb <- de$calc_temp_cor(obs, bad_bad)
  de_bb <- 1 - de$calc_de(obs, bad_bad)
  
  brel_bg <- de$calc_brel_mean(obs, bad_good, sort=sort)
  barea_bg <- de$calc_bias_area(de$calc_brel_res(obs, bad_good, sort=sort))
  r_bg <- de$calc_temp_cor(obs, bad_good)
  de_bg <- 1 - de$calc_de(obs, bad_good)
  
  DE_BB <- list("DE" = de_bb,
                "r" = r_bb,
                "brel" = brel_bb,
                "barea" = barea_bb)
  DE_BG <- list("DE" = de_bg,
                "r" = r_bg,
                "brel" = brel_bg,
                "barea" = barea_bg)
  
  # Create standardized vectors for storing results --------------------------------------
  
  # KGE
  KGE_BB_table <- c("Score" = KGE_BB$KGE.value, 
                    "beta" = KGE_BB$KGE.elements[["Beta"]],
                    "r" = KGE_BB$KGE.elements[["r"]],
                    "alpha" = KGE_BB$KGE.elements[["Alpha"]],
                    "alpha_r" = NA,
                    "r_on_alpha" = NA,
                    "2alpha_r" = NA,
                    "barea" = NA) 
  
  KGE_BG_table <- c("Score" = KGE_BG$KGE.value, 
                    "beta" = KGE_BG$KGE.elements[["Beta"]],
                    "r" = KGE_BG$KGE.elements[["r"]],
                    "alpha" = KGE_BG$KGE.elements[["Alpha"]],
                    "alpha_r" = NA,
                    "r_on_alpha" = NA,
                    "2alpha_r" = NA,
                    "barea" = NA) 
  
  # KGE'
  KGE_m_BB_table <- c("Score" = KGE_m_BB$KGE.value, 
                      "beta" = KGE_m_BB$KGE.elements[["Beta"]],
                      "r" = KGE_m_BB$KGE.elements[["r"]],
                      "alpha" = KGE_m_BB$KGE.elements[["Gamma"]],
                      "alpha_r" = NA,
                      "r_on_alpha" = NA,
                      "2alpha_r" = NA,
                      "barea" = NA) 
  
  KGE_m_BG_table <- c("Score" = KGE_m_BG$KGE.value, 
                      "beta" = KGE_m_BG$KGE.elements[["Beta"]],
                      "r" = KGE_m_BG$KGE.elements[["r"]],
                      "alpha" = KGE_m_BG$KGE.elements[["Gamma"]],
                      "alpha_r" = NA,
                      "r_on_alpha" = NA,
                      "2alpha_r" = NA,
                      "barea" = NA) 
  
  # KGE''
  KGE_m2_BB_table <- c("Score" = KGE_m2_BB$KGE_m2, 
                       "beta" = KGE_m2_BB$Beta_m,
                       "r" = KGE_m2_BB$r,
                       "alpha" = KGE_m2_BB$Alpha,
                       "alpha_r" = NA,
                       "r_on_alpha" = NA,
                       "2alpha_r" = NA,
                       "barea" = NA) 
  
  KGE_m2_BG_table <- c("Score" = KGE_m2_BG$KGE_m2, 
                       "beta" = KGE_m2_BG$Beta_m,
                       "r" = KGE_m2_BG$r,
                       "alpha" = KGE_m2_BG$Alpha,
                       "alpha_r" = NA,
                       "r_on_alpha" = NA,
                       "2alpha_r" = NA,
                       "barea" = NA) 
  
  # non-parametric KGE
  KGENP_BB_table <- c("Score" = KGENP_BB$KGENP, 
                      "beta" = KGENP_BB$Beta,
                      "r" = KGENP_BB$r,
                      "alpha" = KGENP_BB$Alpha,
                      "alpha_r" = NA,
                      "r_on_alpha" = NA,
                      "2alpha_r" = NA,
                      "barea" = NA) 
  
  KGENP_BG_table <- c("Score" = KGENP_BG$KGENP, 
                      "beta" = KGENP_BG$Beta,
                      "r" = KGENP_BG$r,
                      "alpha" = KGENP_BG$Alpha,
                      "alpha_r" = NA,
                      "r_on_alpha" = NA,
                      "2alpha_r" = NA,
                      "barea" = NA) 
  
  ## NSE
  NSE_BB_table <- c("Score" = NSE_BB$NSE, 
                    "beta" = NSE_BB$Beta_m,
                    "r" = NA,
                    "alpha" = NSE_BB$Alpha,
                    "alpha_r" = NA,
                    "r_on_alpha" = NA,
                    "2alpha_r" = NSE_BB$`2_alpha_r`,
                    "barea" = NA) 
  
  NSE_BG_table <- c("Score" = NSE_BG$NSE, 
                    "beta" = NSE_BG$Beta_m,
                    "r" = NA,
                    "alpha" = NSE_BG$Alpha,
                    "alpha_r" = NA,
                    "r_on_alpha" = NA,
                    "2alpha_r" = NSE_BG$`2_alpha_r`,
                    "barea" = NA) 
  
  ## dr
  dr_BB_table <- c("Score" = dr_BB, 
                   "beta" = NA,
                   "r" = NA,
                   "alpha" = NA,
                   "alpha_r" = NA,
                   "r_on_alpha" = NA,
                   "2alpha_r" = NA,
                   "barea" = NA) 
  
  dr_BG_table <- c("Score" = dr_BG, 
                   "beta" = NA,
                   "r" = NA,
                   "alpha" = NA,
                   "alpha_r" = NA,
                   "r_on_alpha" = NA,
                   "2alpha_r" = NA,
                   "barea" = NA) 
  
  ## LME
  LME_BB_table <- c("Score" = LME_BB$LME, 
                    "beta" = LME_BB$Beta,
                    "r" = NA,
                    "alpha" = NA,
                    "alpha_r" = LME_BB$r_alpha,
                    "r_on_alpha" = NA,
                    "2alpha_r" = NA,
                    "barea" = NA) 
  
  LME_BG_table <- c("Score" = LME_BG$LME, 
                    "beta" = LME_BG$Beta,
                    "r" = NA,
                    "alpha" = NA,
                    "alpha_r" = LME_BG$r_alpha,
                    "r_on_alpha" = NA,
                    "2alpha_r" = NA,
                    "barea" = NA) 
  
  ## LCE
  LCE_BB_table <- c("Score" = LCE_BB$LCE, 
                    "beta" = LCE_BB$Beta,
                    "r" = NA,
                    "alpha" = NA,
                    "alpha_r" = LCE_BB$r_alpha,
                    "r_on_alpha" = LCE_BB$r_on_alpha,
                    "2alpha_r" = NA,
                    "barea" = NA) 
  
  LCE_BG_table <- c("Score" = LCE_BG$LCE, 
                    "beta" = LCE_BG$Beta,
                    "r" = NA,
                    "alpha" = NA,
                    "alpha_r" = LCE_BG$r_alpha,
                    "r_on_alpha" = LCE_BG$r_on_alpha,
                    "2alpha_r" = NA,
                    "barea" = NA) 
  
  ## DE
  DE_BB_table <- c("Score" = DE_BB$DE, 
                   "beta" = DE_BB$brel,
                   "r" = DE_BB$r,
                   "alpha" = NA,
                   "alpha_r" = NA,
                   "r_on_alpha" = NA,
                   "2alpha_r" = NA,
                   "barea" = DE_BB$barea) 
  
  DE_BG_table <- c("Score" = DE_BG$DE, 
                   "beta" = DE_BG$brel,
                   "r" = DE_BG$r,
                   "alpha" = NA,
                   "alpha_r" = NA,
                   "r_on_alpha" = NA,
                   "2alpha_r" = NA,
                   "barea" = DE_BG$barea) 
  
  # Create dataframe -------------------------------------------------------------
  
  KGE_results <- data.frame(Param = c("Score", 
                                      "beta",
                                      "r",
                                      "alpha",
                                      "alpha_r",
                                      "r_on_alpha",
                                      "2alpha_r",
                                      "barea"),
                            KGE_BB = KGE_BB_table,
                            KGE_BG = KGE_BG_table,
                            KGE_m_BB = KGE_m_BB_table,
                            KGE_m_BG = KGE_m_BG_table,
                            KGE_m2_BB = KGE_m2_BB_table,
                            KGE_m2_BG = KGE_m2_BG_table,
                            KGENP_BB = KGENP_BB_table,
                            KGENP_BG = KGENP_BG_table,
                            DE_BB = DE_BB_table,
                            DE_BG = DE_BG_table,
                            LME_BB = LME_BB_table,
                            LME_BG = LME_BG_table,
                            LCE_BB = LCE_BB_table,
                            LCE_BG = LCE_BG_table,
                            NSE_BB = NSE_BB_table,
                            NSE_BG = NSE_BG_table,
                            dr_BB = dr_BB_table,
                            dr_BG = dr_BG_table)
  
  # Plot ----------------------------------------------------------------------
  
  # Score plot
  score <- KGE_results |> 
    pivot_longer(ends_with(c("_BB", "_BG"))) |> 
    mutate(criterion = gsub(".{3}$", "", name),
           name = sub(".*(?=.{2}$)", "", name, perl = TRUE)) |> 
    dplyr::filter(Param == "Score") |> 
    select(criterion, name, value) |> 
    mutate(criterion = factor(criterion, levels = rev(c("KGE", "KGE_m", "KGE_m2", "KGENP", "DE",
                                                        "LME", "LCE", "NSE", "dr")))) |> 
    group_by(criterion) |> 
    mutate(minvalue = min(value, na.rm = TRUE))
  
  score_plot <- ggplot(score) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "darkgrey") +
    geom_segment(aes(x = criterion, xend = criterion, y = minvalue, yend = value), 
                 color = "black") +
    geom_point(aes(criterion, value, fill = name), size = 4, shape = 21) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) + 
    scale_x_discrete(labels = c("KGE" = "KGE", 
                                "KGE_m" = expression(paste("KGE","'")), 
                                "KGE_m2" = expression(paste("KGE","''")), 
                                "KGENP" = expression(KGE[NP]), 
                                "DE" = expression(paste("DE","'")),  
                                "LME" = "LME", 
                                "LCE" = "LCE",
                                "NSE" = "NSE", 
                                "dr" = expression(d[r]))) +
    coord_flip(ylim = c(0,1)) +
    scale_fill_manual(name = "",
                      values = c("BB" = "black",
                                 "BG" = "white")) +
    xlab("") +
    ylab("Criterion score") +
    theme_bw() +
    theme(legend.position = "none",
          axis.title.y = element_blank(),
          panel.grid = element_blank(),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.text = element_text(size = 10))
  
  # Parameters plot
  param <- KGE_results |> 
    pivot_longer(ends_with(c("_BB", "_BG"))) |> 
    mutate(criterion = gsub(".{3}$", "", name),
           name = sub(".*(?=.{2}$)", "", name, perl = TRUE)) |> 
    dplyr::filter(Param != "Score") |> 
    mutate(param = case_when(criterion == "KGE" & Param == "beta" ~ "beta",
                             criterion == "KGE" & Param == "alpha" ~ "alpha",
                             criterion == "KGE" & Param == "r" ~ "r",
                             criterion == "KGE_m" & Param == "alpha" ~ "gamma",
                             criterion == "KGE_m2" & Param == "beta" ~ "beta_n",
                             criterion == "KGENP" & Param == "alpha" ~ "alpha_np",
                             criterion == "KGENP" & Param == "r" ~ "r_s",
                             criterion == "DE" & Param == "beta" ~ "brel",
                             criterion == "DE" & Param == "barea" ~ "barea")) |> 
    drop_na() |> 
    select(param, name, value) |> 
    mutate(param = factor(param, levels = rev(c("barea", 
                                                "beta_n", "brel", "beta", 
                                                "alpha", "gamma", "alpha_np", 
                                                "r", "r_s")))) |> 
    group_by(param) |> 
    mutate(minvalue = min(value, na.rm = TRUE))
  
  param_plot <- ggplot(param) +
    geom_segment(x = 10, xend = 6.4, 
                 y = 0, yend = 0, 
                 linetype = "dashed", color = "darkgrey") +
    geom_segment(x = 6.6, xend = 0, 
                 y = 1, yend = 1, 
                 linetype = "dashed", color = "darkgrey") +
    geom_segment(aes(x = param, xend = param, y = minvalue, yend = value), 
                 color = "black") +
    geom_point(aes(param, value, fill = name), size = 4, shape = 21) +
    scale_y_continuous(breaks = c(-0.5, 0, 0.5, 1)) +
    scale_x_discrete(labels = c("beta" = expression(beta), 
                                "beta_n" = expression(beta[n]), 
                                "brel" = expression(B[rel]), 
                                "alpha" = expression(alpha), 
                                "gamma" = expression(gamma),  
                                "alpha_np" = expression(alpha[NP]), 
                                "r" = expression(r),
                                "r_s" = expression(r[s]), 
                                "barea" = expression(B[area]))) +
    coord_flip() +
    scale_fill_manual(name = "",
                      values = c("BB" = "black",
                                 "BG" = "white")) +
    xlab("") +
    ylab("Parameter value") +
    theme_bw() +
    theme(legend.position = "none",
          axis.title.y = element_blank(),
          panel.grid = element_blank(),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.text = element_text(size = 10))
  
  # Get legend
  legend <- get_legend(
    ggplot(param) +
      geom_point(aes(param, value, fill = name), size = 4, shape = 21) +
      scale_fill_manual(name = "", 
                        values = c("BB" = "black", "BG" = "white"),
                        labels = c("BB" = "BB model",
                                   "BG" = "BG model")) +
      theme_bw() +
      theme(legend.position = "bottom",
            axis.title = element_text(size = 12),
            axis.text = element_text(size = 10),
            legend.text = element_text(size = 10)))
  
  # Final plot
  plot_grid(
    plot_grid(score_plot, param_plot),
    legend,
    ncol = 1, 
    rel_heights = c(1, 0.1)
  )
  
}