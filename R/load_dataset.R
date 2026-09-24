load_ds <- function(x) {
  fname <- basename(tools::file_path_sans_ext(x))
  fname <- gsub("([0-9]{4})", " (\\1)", fname)
  out <- mget(load(x, envir = (tmp_env <- new.env())), envir = tmp_env)[[1]]

  if (!is.null(out$controls) && !length(out$controls) == 0L) {
    out$X <- model.matrix(~ ., data = out$controls)
  }

  out$name <-fname 
  x <- out 

  out  <- tibble(
    name = x$name, 
    iv_diag = list(x$iv_diag), 
    f_effective = x$iv_diag$F_stat["F.effective"],
    w = list(x$w),
    w_name = x$w_name,
    z = list(x$z),
    z_name = x$z_name,
    y = list(x$y),
    y_name = x$y_name,
    controls = list(x$controls), 
    controls_name = list(x$controls_name), 
    weights = list(x$weights),
    weights_name = x$weights_name,
    FE = list(x$FE),
    FE_name = list(x$FE_name),
    X = list(x$X), 
    df = list(x$df), 
  )

  stopifnot(nrow(out) == 1)

  return(out)
}

