make_contour_plot <- function(x, rvs, bench, grid_resolution = 400) {
  if (!is.null(x$X)) {
      vars <- cbind(x$w, x$z, x$y)
      residuals <- stats::lm.fit(
          x = x$X,
          y = vars,
          singular.ok = TRUE
        )$residuals

      x$w <- residuals[,1]
      x$z <- residuals[,2]
      x$y <- residuals[,3]
  }

  soo_rv = round(rvs$rv[1],3)
  iv_rv = round(rvs$rv[2],3)

  cmat <- ivrobustness:::CovarianceMatrix$new(
    x$w, 
    x$z, 
    x$y
  )

  extended_cmat <- cmat$extend(c(0,0,0))

  tau <-extended_cmat$tau() 
  iv_est <-extended_cmat$iv_estimate() 
  soo_est <-extended_cmat$soo_estimate() 

  get_new_tau <- function(x,y) { 
    extended_cmat <- cmat$extend(c(0,x,y))
    extended_cmat$tau()
  }

  df <- expand_grid(
    r_1 = seq(-.99,.99, length.out = grid_resolution),
    r_2 = seq(-.99,.99, length.out = grid_resolution)
  ) %>% 
    mutate(
      true_te = map2(r_1, r_2, safely(get_new_tau)),
      result = map_dbl(true_te, "result"),
      soo_bias = soo_est - result,
      iv_bias = iv_est - result,
      error = map(tau, "error"),
      true_te = NULL
    )  

  out_plot <- df %>%
    mutate(iv_advantage = abs(soo_bias) - abs(iv_bias)) %>%
    ggplot(aes(x = r_1, y = r_2)) +
    geom_raster(aes(fill = iv_advantage)) +
    geom_contour(
      aes(z = iv_advantage,
          colour = "Equal absolute bias",
          linetype = "Equal absolute bias"),
      breaks = 0, linewidth = 0.8, na.rm = TRUE
    ) +
    geom_contour(
      aes(z = iv_bias,
          colour = "Zero IV bias",
          linetype = "Zero IV bias"),
      breaks = 0, linewidth = 0.9, na.rm = TRUE
    ) 

    if (length(bench) != 0L)  {
      out_plot  <- out_plot + 
      geom_point(
        data = bench,
        aes(x = rho_2, y = rho_3),
        shape = 21, size = 2.8, fill = "white", colour = "grey15"
      ) +
      geom_text_repel(
        data = bench,
        aes(x = rho_2, y = rho_3, label = variable),
        size = 3.4, seed = 1, 
        box.padding = 0.4, min.segment.length = 0,
        colour = "grey15", segment.color = "grey50"
      )
    }


    out_plot <- out_plot + scale_fill_gradient2(
      name = "|SOO bias| - |IV bias|",
      low = "#E8B89B", mid = "white", high = "#99C5DF",
      midpoint = 0, na.value = "grey90"
    ) +
    scale_colour_manual(
      name = NULL,
      values = c("Equal absolute bias" = "grey15",
                "Zero IV bias" = "#6A3D9A")
    ) +
    scale_linetype_manual(
      name = NULL,
      values = c("Equal absolute bias" = "solid",
                "Zero IV bias" = "dashed")
    ) +
    coord_equal() +
    labs(
      title = str_glue("IV versus selection on observables {x$name}"),
      subtitle = str_glue("SOO RV: {soo_rv}, IV RV: {iv_rv}"),
      caption = "Blue: IV less biased · Orange: SOO less biased",
      x = expression(rho[UZ * "|" * W]),
      y = expression(rho[UY * "|" * W * "," * Z])
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold"),
      legend.position = "bottom",
      legend.box = "vertical"
    )

  return(list(out_plot))
}


save_contour_plot <- function(x, plt) {
  out <- str_glue("figures/contours/{x$name}.png")

  ggsave(plot = plt, out, width = 7, height = 7)

  return(out)
}
