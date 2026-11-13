# _targets.R file
library(targets)
library(crew)
library(tarchetypes)

tar_source()

tar_option_set(
  packages = c("tidyverse", "ivrobustness", "fixest", "ggrepel", "MASS", "patchwork", "ggcorrplot"),
  controller = crew_controller_local(workers = 12)
)

make_bens_bound_plot <- function(df, rvs) {
  rvs <- rvs %>%
    filter(method == "Iv")

  out <- df %>% 
    mutate(
      rf_rs = pmap_dbl(list(w, y, X), function(w, y, X) {
        if (!is.null(X)) {
          controls <- cbind(1, as.matrix(X))
          res <- lm.fit(x = controls, y = cbind(w, y))$residuals

          return(cor(res[, 1], res[, 2])^2)
        } else {
          return(cor(w,y)^2)
        }
      })
    )  %>%
    inner_join(rvs)

  median(out$rv / out$rf_rs)


  out %>% 
    ggplot(aes(x=rf_rs, y=rv)) + 
    geom_point() + 
    geom_abline(slope = 1, linetype = 2) + 
    coord_equal() + 
    ggtitle(
      "Robustness Value vs Reduced Form R-Squared",
      subtitle = "69 IV Analyses from Lal et al (2024)"
    ) + 
    xlab("Reduced Form R-sqared") + 
    ylab("IV Robustness Value")

  out_path = "RF_bound.svg"
  ggsave(out_path,  width = 6, height = 4)
  ggsave("RF_bound.png", width = 6, height = 4)

  return(out_path)
}



calc_rvs <- function(df) { 
  stopifnot(nrow(df) == 1)
  out_df <- robustness_values(df$w[[1]], df$z[[1]], df$y[[1]], df$X[[1]])

  out_df$name <- df$name

  return(out_df)
}

calc_benches <- function(df) { 
  stopifnot(nrow(df) == 1)

  w <- as.double(haven::zap_labels(df$w[[1]]))
  z <- as.double(haven::zap_labels(df$z[[1]]))
  y <- as.double(haven::zap_labels(df$y[[1]]))
  X <- df$X[[1]]

  if (is.null(X) || length(X) == 0L) {
    return(tibble())
  }

  out_df <- benchmark_covariates(
    w, z, y, X
  )


  out_df$name <- df$name

  return(out_df)
}

list(
  tar_target(files, list.files( "datasets", full.names = T)),
  tar_target(datasets, subtract_fes(load_ds(files)), pattern = map(files)), 

  tar_target(summary_table, make_summary_table(datasets)), 
  tar_target(exaggeration_plot, make_exaggeration_plot(summary_table)), 


  tar_target(sim_df, expand_grid(
    r_yz = seq(-.95,.95, length.out = 20),
    r_zw = seq(-.95,.95, length.out = 20),
    p_yw_z = seq(-.95,.95, length.out = 20),
  )),

  tar_target(sims, 
    run_one_sim(sim_df), 
    pattern = map(sim_df), 
    error = "null"
  ),

  tar_target(rvs, calc_rvs(datasets), pattern = map(datasets)), 
  tar_target(benches, calc_benches(datasets), pattern = map(datasets)), 
  tar_target(analysis_df, full_join(datasets, rvs, by = "name") %>% full_join(benches, by = "name")),

  tar_target(scatterplots, make_scatterplots(datasets, rvs), format = "file"),
  tar_target(sim_scatterplots, make_sim_scatters(sims), format = "file"),
  tar_target(bens_bound_plot, make_bens_bound_plot(datasets, rvs), format = "file"),

  tar_target(iv_contour_plots, make_iv_contour_plot(datasets), pattern = map(datasets), format = "file"),
  tar_target(contour_plots, make_contour_plot(datasets, rvs, benches), pattern = map(datasets, rvs, benches)),
  tar_target(saved_contour_plots, save_contour_plot(datasets, contour_plots), pattern = map(datasets, contour_plots)),
  tar_quarto(report, "report2.qmd")
)

