library(tidyverse)
library(targets)

datasets <- targets::tar_read(datasets)
df <- datasets[2,]

calc_rvs(df)



process_ds <- function(df) {
  stopifnot(nrow(df) == 1)

  x <- as.list(df)
  x <- lapply(x, `[[`, 1L)

  if (!is.null(x$FE) && length(x$FE) > 0L) {
    vars <- cbind(w = x$w, z = x$z, y = x$y, x$X)

    vars <- fixest::demean(
      X = vars,
      f = x$FE,
      na.rm = FALSE,
      as.matrix = TRUE
    )
  } 


  df$w <- list(as.double(vars[, 1]))
  df$z <- list(as.double(vars[, 2]))
  df$y <- list(as.double(vars[, 3]))
  df$X <- list(vars[, -seq_len(3), drop = FALSE])

  df
}

process_ds(datasets[11,])

