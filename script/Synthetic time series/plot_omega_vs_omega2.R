#' Plot omega vs omega2 with defined sampling of omega parameter
#'
#' @param data Numeric vector of discharge time series
#' @param perf_c String vector of performance criteria. Can be kge, kge_m, kgenp, nse, kge_m2, de, lme
#' @param omega_log_min Numeric value. Minimum threshold for omega sampling
#' @param omega_log_max Numeric value. Maximum threshold for omega sampling
#' @param step Numeric value. Step of the sampling between min and max
#'
#' @return A plot
#' @export
#'
#' @examples

plot_omega_vs_omega2 <- function(data,
                                 perf_c = c("kge", "kge_m", "kgenp", "nse"),
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
    
    # Get max values
    for (crit in perf_c) {
      all_c[[paste0(i)]][[paste0(crit, "_max_omega2")]] <- 
        max(dplyr::filter(all_c[[paste0(i)]]$param, get(crit) == max(get(crit)))$omega)
    }
  }
  
  # Get final dataset
  syn_df <- data.frame()
  # because if no filter omega returns numeric(0)
  name_filter <- paste0(perf_c, "_max_omega2") 
  
  for (i in seq(length(all_c))) {
    clear_list <- all_c[[i]][c("omega", name_filter)]
    clear_list_df <- as.data.frame(clear_list)
    syn_df <- rbind(syn_df, clear_list_df)
  }
  
  # Plot results
  plot_values <- c("kge_max_omega2" = cblind_bw_palette[2],
                   "kge_m_max_omega2" = cblind_bw_palette[3],
                   "kge_m2_max_omega2" = cblind_bw_palette[5],
                   "kgenp_max_omega2" = cblind_bw_palette[4],
                   "de_max_omega2" = cblind_bw_palette[6],
                   "lme_max_omega2" = cblind_bw_palette[7],
                   "lce_max_omega2" = cblind_bw_palette[8],
                   "nse_max_omega2" = cblind_bw_palette[1],
                   "dr_max_omega2" = cblind_bw_palette[1])
  
  plot_labels <- c("kge_max_omega2" = expression(KGE),
                   "kge_m_max_omega2" = expression(paste("KGE","'")),
                   "kge_m2_max_omega2" = expression(paste("KGE","''")),
                   "kgenp_max_omega2" = expression(KGE[NP]),
                   "de_max_omega2" = expression(paste("DE","'")),
                   "lme_max_omega2" = expression(LME),
                   "lce_max_omega2" = expression(LCE),
                   "nse_max_omega2" = expression(NSE),
                   "dr_max_omega2" = expression(d[r]))
  
  long_syn_df <- syn_df |> 
    pivot_longer(name_filter, names_to = "criterion", values_to = "omega2")
  
  long_syn_df |> 
    ggplot(aes(omega, omega2, color = criterion)) +
    geom_line(alpha = 0.8) +
    coord_cartesian(xlim = c(min(long_syn_df$omega), max(long_syn_df$omega)), 
                    ylim = c(min(long_syn_df$omega), max(long_syn_df$omega))) +
    scale_x_continuous(trans = "log10", 
                       expand = expansion(mult = c(0, 0)),
                       breaks = c(0.5, 1, 2)) +
    scale_y_continuous(trans = "log10", 
                       breaks = c(0.5, 1, 2)) +
    scale_color_manual(name = "Best model", 
                       values = plot_values[names(plot_values) %in% name_filter],
                       labels = plot_labels[names(plot_labels) %in% name_filter]) +
    xlab(expression(omega[1])) +
    ylab(expression(omega[2])) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          legend.text = element_text(size = 10, hjust = 0),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10))
}



