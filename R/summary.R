make_summary_table <- function(datasets) { 
  iv_diag_to_tibble <- function(x) {
    tibble(
      ols_coefficient = x$est_ols[1,1],
      iv_coefficient = x$est_2sls[1,1],
      rf_coefficient = x$est_rf[,'Coef'],
      fs_coefficient = x$est_fs[,'Coef'],
      N = x$N
    )
  }

  datasets %>% 
    mutate(
      iv_df = map(iv_diag, iv_diag_to_tibble)
    ) %>% 
    dplyr::select(name, f_effective, iv_df) %>% 
    unnest(iv_df)
}

make_exaggeration_plot <- function(summary_table) { 
  summary_table %>% 
    ggplot(aes(x = abs(fs_coefficient), y = abs(iv_coefficient))) + 
    geom_point() +  
    scale_x_log10() + 
    scale_y_log10() + 
    geom_smooth(method = "lm") +
    xlab("Absolute Value of First-Stage Coefficient") + 
    ylab("Absolute Value of IV Coefficient")  + 
    labs(
      caption = "Negative Relationship Between FS and 2SLS coefficients"
    ) + 
    theme_minimal()

   out_path <- "figures/summary_figures/exaggeration.svg"
   ggsave(out_path)

   return(out_path)
}

