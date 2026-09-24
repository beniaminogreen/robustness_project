subtract_fes <- function(df, standardize = TRUE) {
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

    df$w <- list(as.double(vars[, 1]))
    df$z <- list(as.double(vars[, 2]))
    df$y <- list(as.double(vars[, 3]))
    df$X <- list(vars[, -seq_len(3), drop = FALSE])
  } 

  if (standardize) { 
    df$y[[1]]  <- df$y[[1]] / sd(df$y[[1]])
    df$z[[1]]  <- df$z[[1]] / sd(df$z[[1]])
    df$w[[1]]  <- df$w[[1]] / sd(df$w[[1]])
  }

  df
}
