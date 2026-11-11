make_sigma <- function(r12, r13, p23_1) {
  stopifnot(all(abs(c(r12, r13, p23_1)) < 1))

  r23 <- r12 * r13 +
    p23_1 * sqrt((1 - r12^2) * (1 - r13^2))

  matrix(
    c(1,   r12, r13,
      r12, 1,   r23,
      r13, r23, 1),
    nrow = 3,
    byrow = TRUE
  )
}

make_dataset <- function(r_yz, r_zw, p_yw_z, n = 10^6) {
  sigma <- make_sigma(r_yz, r_zw, p_yw_z)

  X <- matrix(rnorm(n * 3), nrow = n) %*% chol(sigma)
  colnames(X) <- c("Z", "Y", "W")

  as_tibble(X)
}

run_one_sim <- function(config) { 
  stopifnot(nrow(config) == 1)
  df <- make_dataset(config$r_yz[1], config$r_zw[1], config$p_yw_z[1])
  rvs <- robustness_values(df$W, df$Z, df$Y)
  return(bind_cols(config, rvs))
}

