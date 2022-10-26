#' Plot synthetic model distribution with defined sampling of omega parameter
#'
#' @param data Numeric vector of discharge time series
#' @param perf_c String vector of performance criteria. Can be kge, kge_m, kgenp, nse, kge_m2, de, lme, dr
#' @param ncol Numeric value. Number of column in the plot
#' @param omega_log_min Numeric value. Minimum threshold for omega sampling
#' @param omega_log_max Numeric value. Maximum threshold for omega sampling
#' @param step Numeric value. Step of the sampling between min and max values
#'
#' @return A plot
#' @export
#'
#' @examples

plot_model_distribution <- function(data,
                                    perf_c = c("kge", "kge_m", "kgenp", "nse"),
                                    ncol = 2,
                                    omega_log_min = -0.36,
                                    omega_log_max = 0.36,
                                    step = 0.01) {
  
  # Isolate syn1
  syn1_q <- syn1$discharge
  syn2_q <- c(syn1_q, syn1_q)
  
  # Create omega  parameters
  omega <- seq(omega_log_min, omega_log_max, step)
  omega <- 10^omega
  
  # All combination of unique parameters
  param <- data.frame(omega = omega) |> 
    mutate(kge = NA,
           kge_m = NA,
           kgenp = NA,
           nse = NA,
           dr = NA,
           kge_m2 = NA,
           lme = NA,
           lce = NA,
           de = NA)
  
  # Create all translations
  all_t <- list()
  
  for (i in seq(nrow(param))) {
    all_t[[i]] <- syn1_q * param$omega[i]
  }
  
  # Create all combinations
  all_c <- list()
  
  for (i in seq(length(all_t))) {
    print(i)
    all_c[[paste0(i)]][["omega"]] <- param$omega[i]
    all_c[[paste0(i)]][["param"]] <- param
    
    for (j in seq(length(all_t))) {
      all_c[[paste0(i)]][[paste0("x", j)]] <- c(all_t[[i]], all_t[[j]])
      
      # KGE
      if ("kge" %in% perf_c) {
        all_c[[paste0(i)]]$param$kge[j] <- 
          KGE(all_c[[paste0(i)]][[paste0("x", j)]], syn2_q, method = "2009")
      }
      
      # KGE'
      if ("kge_m" %in% perf_c) {
        all_c[[paste0(i)]]$param$kge_m[j] <- 
          KGE(all_c[[paste0(i)]][[paste0("x", j)]], syn2_q, method = "2012")
      }
      
      # KGENP
      if ("kgenp" %in% perf_c) {
        all_c[[paste0(i)]]$param$kgenp[j] <- 
          (RNP(all_c[[paste0(i)]][[paste0("x", j)]], syn2_q))$KGENP
      }
      
      # NSE
      if ("nse" %in% perf_c) {
        all_c[[paste0(i)]]$param$nse[j] <- 
          NSE(all_c[[paste0(i)]][[paste0("x", j)]], syn2_q)
        
      }
      
      # dr
      if ("dr" %in% perf_c) {
        all_c[[paste0(i)]]$param$dr[j] <- 
          hydroErr$dr(all_c[[paste0(i)]][[paste0("x", j)]], syn2_q)
      }
      
      # KGE''
      if ("kge_m2" %in% perf_c) {
        all_c[[paste0(i)]]$param$kge_m2[j] <- 
          KGE_m2(all_c[[paste0(i)]][[paste0("x", j)]], syn2_q)
      }
      
      # DE
      if ("de" %in% perf_c) {
        all_c[[paste0(i)]]$param$de[j] <- 
          1 - de$calc_de(syn2_q, all_c[[paste0(i)]][[paste0("x", j)]])
      }
      
      # LME
      if ("lme" %in% perf_c) {
        all_c[[paste0(i)]]$param$lme[j] <- 
          LME(all_c[[paste0(i)]][[paste0("x", j)]], syn2_q)
        
      }
      
      # LCE
      if ("lce" %in% perf_c) {
        all_c[[paste0(i)]]$param$lce[j] <- 
          LCE(all_c[[paste0(i)]][[paste0("x", j)]], syn2_q)
      }
    } 
    
    # Get max, min, good values
    for (crit in perf_c) {
      all_c[[paste0(i)]][[paste0("max_", crit)]] <- 
        max(dplyr::filter(all_c[[paste0(i)]]$param, omega != 1)[[crit]])
      all_c[[paste0(i)]][[paste0("min_", crit)]] <- 
        min(all_c[[paste0(i)]]$param[[crit]])
      all_c[[paste0(i)]][[paste0("good_", crit)]] <- 
        dplyr::filter(all_c[[paste0(i)]]$param, omega == 1)[[crit]]
    }
  }
  
  # Get final dataset
  syn_df <- data.frame()
  
  for (i in seq(length(all_c))) {
    clear_list <- all_c[[i]][c("omega", 
                               "max_kge", "min_kge", "good_kge", 
                               "max_kge_m", "min_kge_m", "good_kge_m", 
                               "max_kgenp", "min_kgenp", "good_kgenp", 
                               "max_nse", "min_nse", "good_nse", 
                               "max_dr", "min_dr", "good_dr", 
                               "max_kge_m2", "min_kge_m2", "good_kge_m2", 
                               "max_de", "min_de", "good_de", 
                               "max_lme", "min_lme", "good_lme", 
                               "max_lce", "min_lce", "good_lce")]
    clear_list_df <- as.data.frame(clear_list)
    syn_df <- rbind(syn_df, clear_list_df)
  }
  
  # Generate plot
  
  ## Arrange max
  max_syn_df <- syn_df |> 
    select(omega, 
           max_kge, max_kge_m, max_kgenp, 
           max_nse, max_dr, max_kge_m2, 
           max_de, max_lme, max_lce) |> 
    pivot_longer(c(starts_with("max")), names_to = "group", values_to = "max") |> 
    mutate(group = case_when(group == "max_kge" ~ "kge",
                             group == "max_kge_m" ~ "kge_m",
                             group == "max_kgenp" ~ "kgenp",
                             group == "max_nse" ~ "nse",
                             group == "max_dr" ~ "dr",
                             group == "max_kge_m2" ~ "kge_m2",
                             group == "max_de" ~ "de",
                             group == "max_lme" ~ "lme",
                             group == "max_lce" ~ "lce"))
  
  ## Arrange min
  min_syn_df <- syn_df |> 
    select(omega, 
           min_kge,  min_kge_m, min_kgenp, 
           min_nse, min_dr, min_kge_m2, 
           min_de, min_lme, min_lce) |> 
    pivot_longer(c(starts_with("min")), names_to = "group", values_to = "min") |> 
    mutate(group = case_when(group == "min_kge" ~ "kge",
                             group == "min_kge_m" ~ "kge_m",
                             group == "min_kgenp" ~ "kgenp",
                             group == "min_nse" ~ "nse",
                             group == "min_dr" ~ "dr",
                             group == "min_kge_m2" ~ "kge_m2",
                             group == "min_de" ~ "de",
                             group == "min_lme" ~ "lme",
                             group == "min_lce" ~ "lce"))
  
  ## Arrange good
  good_syn_df <- syn_df |> 
    select(omega, 
           good_kge,  good_kge_m, good_kgenp, 
           good_nse, good_dr, good_kge_m2, 
           good_de, good_lme, good_lce) |> 
    pivot_longer(c(starts_with("good")), names_to = "group", values_to = "good") |> 
    mutate(group = case_when(group == "good_kge" ~ "kge",
                             group == "good_kge_m" ~ "kge_m",
                             group == "good_kgenp" ~ "kgenp",
                             group == "good_nse" ~ "nse",
                             group == "good_dr" ~ "dr",
                             group == "good_kge_m2" ~ "kge_m2",
                             group == "good_de" ~ "de",
                             group == "good_lme" ~ "lme",
                             group == "good_lce" ~ "lce"))
  
  # Combine max, min, good
  full_syn_df <- good_syn_df |> 
    left_join(max_syn_df, by = c("omega", "group")) |> 
    left_join(min_syn_df, by = c("omega", "group"))
  
  # Calculate performance for syn example and create df for geom point of example
  syn2_model2 <- c(syn1_q * 0.75, syn1_q)
  
  point_syn_example <- data.frame(omega = 0.75,
                                  group = perf_c) |> 
    mutate(good = case_when(group == "kge" ~ KGE(syn2_model2, syn2_q, method = "2009"),
                            group == "kge_m" ~ KGE(syn2_model2, syn2_q, method = "2012"),
                            group == "kgenp" ~ (RNP(syn2_model2, syn2_q))$KGENP,
                            group == "nse" ~ NSE(syn2_model2, syn2_q),
                            group == "dr" ~ hydroErr$dr(syn2_model2, syn2_q),
                            group == "kge_m2" ~ KGE_m2(syn2_model2, syn2_q),
                            group == "de" ~ 1 - de$calc_de(syn2_q, syn2_model2),
                            group == "lme" ~ LME(syn2_model2, syn2_q),
                            group == "lce" ~ LCE(syn2_model2, syn2_q)))
  
  # Final plot
  full_syn_df |> 
    drop_na() |> 
    ggplot(aes(x = omega, y = good)) +
    geom_ribbon(aes(ymin = min, 
                    ymax = max, 
                    group = group, 
                    fill = "Distribution of models with \u03C92≠1"),
                alpha = 0.3,
                color = "black",
                linetype = "dotted",
                size = 0.3) +
    geom_line(aes(color = "Models with \u03C92=1")) +
    geom_point(data = point_syn_example, 
               aes(shape = '"Bad-Good" model shown in Figure 2'), 
               color = "black", 
               fill = "white", 
               size = 2) +
    facet_wrap(~ factor(group, levels = c("kge", "kge_m", "kge_m2", 
                                          "nse", "kgenp", "de",
                                          "dr", "lme", "lce")), 
               ncol = ncol, 
               labeller = as_labeller(c("kge" = "KGE", 
                                        "kge_m" = paste(expression(paste("KGE",
                                                                         "'"))), 
                                        "kgenp" = "KGE[NP]", 
                                        "nse" = "NSE", 
                                        "dr" = "d[r]", 
                                        "kge_m2" = paste(expression(paste("KGE",
                                                                          "''"))), 
                                        "de" = paste(expression(paste("DE",
                                                                      "'"))), 
                                        "lme" = "LME", 
                                        "lce" = "LCE"), label_parsed)) +
    coord_cartesian(xlim = c(min(full_syn_df$omega), max(full_syn_df$omega)), 
                    ylim = c(-0.41, 1)) +
    scale_x_continuous(trans = "log10", 
                       expand = expansion(mult = c(0, 0)),
                       breaks = c(0.5, 1, 2)) +
    scale_fill_manual(name = "", 
                      values = "#ffffd6",
                      labels = expression(paste(Distribution~of~models, 
                                                " for ", 
                                                all~omega[2]~values))) +
    scale_color_manual(name = "", 
                       values = "#696969",
                       labels = expression(paste(Models~with~omega[2], 
                                                 " = ", 
                                                 1))) +
    scale_shape_manual(name = "", values = 21) +
    xlab(expression(omega[1])) +
    ylab("Criterion score") +
    theme_bw() +
    theme(legend.position = "bottom",
          panel.grid.minor = element_blank(),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          strip.text = element_text(size = 10),
          legend.text = element_text(size = 10)) +
    guides(shape = guide_legend(order = 1), 
           color = guide_legend(order = 2), 
           fill = guide_legend(order = 3))
}



