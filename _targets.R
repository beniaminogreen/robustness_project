# _targets.R file
library(targets)
library(crew)
library(tarchetypes)

tar_source()

tar_option_set(
  packages = c("tidyverse", "ivrobustness", "fixest", "ggrepel", "MASS"),
  controller = crew_controller_local(workers = 12)
)


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

  tar_target(sim_df, expand_grid(
    r_yz = seq(.05,.95, length.out = 19),
    r_zw = seq(.05,.95, length.out = 19),
    p_yw_z = seq(.05,.95, length.out = 19),
  )),

  tar_target(sims, 
    run_one_sim(sim_df), 
    pattern = map(sim_df), 
    error = "null"
  ),

  tar_target(rvs, calc_rvs(datasets), pattern = map(datasets)), 
  tar_target(benches, calc_benches(datasets), pattern = map(datasets)), 
  tar_target(analysis_df, full_join(datasets, rvs, by = "name") %>% full_join(benches, by = "name")),

  tar_target(iv_contour_plots, make_iv_contour_plot(datasets), pattern = map(datasets), format = "file"),
  tar_target(contour_plots, make_contour_plot(datasets, rvs, benches), pattern = map(datasets, rvs, benches)),
  tar_target(saved_contour_plots, save_contour_plot(datasets, contour_plots), pattern = map(datasets, contour_plots)),
  tar_quarto(report, "report.qmd")
)

