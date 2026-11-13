make_scatterplots <- function(datasets, rvs) { 
  df <- full_join(datasets, rvs, by = "name")

  experimental_studies <- c(
    "Alt (2015)",
    "Blattman (2014)",
    "Chong (2019)",
    "Coppock (2016)",
    "Delao (2013)",
    "Gerber (2010)",
    "Healy (2013)",
    "Lopez (2020)",
    "McClendon (2014)",
    "Pianzola (2019)",
    "Webster (2022)",
    "Wood (2022)"
  )


  df <- df %>% 
    mutate(
      experimental = name %in% experimental_studies,
      selected_iv_method = ifelse(experimental, "ExogenousIV","Iv"), 
      first_stage = map2_dbl(z,w, function(z,w){cov(z,w)/var(z)}) 
    ) %>% 
    filter(
      method == "Soo" | method == selected_iv_method
    ) %>% 
    mutate(
      method = ifelse(method == "Soo", "SOO", "IV")
    )

  scatter_df <- df %>% 
    dplyr::select(name, method, rv, experimental) %>% 
    pivot_wider(id_cols = c(name, experimental), names_from = method, values_from = rv) 

  scatter_df %>% 
    ggplot(aes(x = SOO, y = IV)) + 
    geom_point(aes(col = experimental)) + 
    geom_text_repel(aes(label = name)) + 
    xlim(0,1) + 
    ylim(0,1) + 
    xlab("SOO Robustness Value") + 
    ylab("IV Robustness Value")  + 
    theme_minimal() + 
    labs(
      col = "Experimetntal Study?"
    ) + 
    ggtitle(
      "IV vs SOO Robustness Values", 
      subtitle = "69 IV Analyses from Lal et. al (2024)"
    )
  
    ggtitle("Scatterplot of Robustness Values")

  bivariate_1 <- "figures/bivariate_1.svg"
  ggsave(bivariate_1, height = 4)

  scatter_df %>% 
    ggplot(aes(x = SOO, y = IV)) + 
    geom_point(aes(col = experimental)) + 
    geom_text_repel(aes(label = name)) + 
    xlab("SOO Robustness Value") + 
    ylab("IV Robustness Value")  + 
    theme_minimal() + 
    labs(
      col = "Experimetntal Study?"
    )  + 
    ggtitle(
      "IV vs SOO Robustness Values", 
      subtitle = "69 IV Analyses from Lal et. al (2024)"
    )

    ggtitle("Scatterplot of Robustness Values")

  bivariate_2 <- "figures/bivariate_2.svg"
  ggsave(bivariate_2,  height = 4)



  df <- df %>% 
    mutate(
      tau = pmap_dbl(list(w, z , y), function(w,z,y) { 
      z_hat <- predict(lm(z ~ w))
      model <- lm(y ~ z_hat)
      return(coef(model)[2])
    }), 
    first_stage = pmap_dbl(list(w, z , y), function(w,z,y) { 
      fs <- coef(lm(z ~ w))[2]
      return(fs)
    }) 
    )

  plot_df <- df %>% 
    filter(method == "IV")


  rsquared_1 <- round(cor(log(abs(plot_df$tau)), log(plot_df$rv))^2,3)
  rsquared_2 <- round(cor(log(abs(plot_df$f_effective)), log(plot_df$rv))^2,3)
  rsquared_3 <- round(cor(log(abs(plot_df$first_stage)), log(plot_df$rv))^2,3)

  plot_df %>% 
    ggplot(aes(x = abs(tau), y = rv)) + 
    geom_point() + 
    geom_smooth(method = "lm") + 
    theme_bw()  + 
    scale_x_log10() + 
    scale_y_log10() + 
    xlab("Absolute Value of Effect Size in Standard Deviations (Log Scale)") + 
    ylab("IV Robustness Value (Log Scale)") + 
    ggtitle(
      "IV Robustness Value vs 2SLS Coefficient",
      subtitle = str_glue("R-Squared : {rsquared_1}")
    ) 
l
  path_1 <- "figures/rv_vs_effect_estimate.svg"
  ggsave(path_1, height = 3)

  plot_df %>% 
    ggplot(aes(x = f_effective, y = rv)) + 
    geom_point() + 
    geom_smooth(method = "lm") + 
    theme_bw()  + 
    scale_x_log10() + 
    scale_y_log10() + 
    xlab("Effective First-Stage F-Statistic (Log Scale)") + 
    ylab("IV Robustness Value (Log Scale)") + 
    ggtitle(
      "IV Robustness Value vs F-Statistic",
      subtitle = str_glue("R-Squared : {rsquared_2}")
    ) 

  path_2 <- "figures/rv_vs_f_stat.svg"
  ggsave(path_2, height = 3)


  plot_df %>% 
    ggplot(aes(x = abs(first_stage), y = rv)) + 
    geom_point() + 
    geom_smooth(method = "lm") + 
    theme_bw()  + 
    scale_x_log10() + 
    scale_y_log10() + 
    xlab("Effective First-Stage F-Statistic (Log Scale)") + 
    ylab("IV Robustness Value (Log Scale)") + 
    ggtitle(
      "IV Robustness Value vs First-Stage Coefficient",
      subtitle = str_glue("R-Squared : {rsquared_3}")
    ) 

  path_3 <- "figures/rv_vs_first_stage.svg"
  ggsave(path_3, height = 3)


  return(c(bivariate_1, bivariate_2,path_1, path_2, path_3))
}

make_sim_scatters <- function(sims) {
  data <- sims %>% mutate(
    r2_wy = (r_yz * r_zw + p_yw_z * sqrt((1 - r_yz^2) * (1 - r_zw^2)))^2
  )

  data %>% 
    filter(method == "Iv")  %>% 
    ggplot(aes(x=r2_wy, y = rv)) +
    geom_abline(slope = 1, linetype = 2) + 
    geom_point() + 
    theme_bw()  + 
    ggtitle(
      "IV-RV by First-Stage R-Squared",
    )  +
    xlab("Reduced-Form R-Squared") + 
    ylab("IV Robustness Value") 
  out_1 <- "figures/sim_scatter_rf_vs_rv.png"
  ggsave(out_1, height = 4)

  # data %>% 
  #   filter(method == "ExogenousIV")  %>% 
  #   ggplot(aes(x=r2_wy, y = rv)) +
  #   geom_abline(slope =1 ) + 
  #   geom_point(alpha =.4) + 
  #   theme_bw()  + 
  #   ggtitle(
  #     "IV-RV by First-Stage R-Squared",
  #     "Imposing Exogeneity"
  #   )  +
  #   xlab("First-Stage R-Squared") + 
  #   ylab("Robustness Value, Imposing Exogeneity of Instrument")

  data %>% 
    filter(method == "ExcludableIV")  %>% 
    ggplot(aes(x=r2_wy, y = rv)) +
    geom_abline(slope =1 ) + 
    geom_point() + 
    theme_bw()  + 
    ggtitle(
      "IV-RV by First-Stage R-Squared",
      "Imposing Exclusion Restriction"
    )   +
    xlab("First-Stage R-Squared") + 
    ylab("IV Robustness Value,\nImposing Exclusion Restriction")
  out_2 <- "figures/sim_scatter_rf_vs_rv_er.png"
  ggsave(out_2, height = 4)


  return(c(out_1, out_2))
}
