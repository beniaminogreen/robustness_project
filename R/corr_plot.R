make_corr_plot <- function(W, Z, Y) { 
  ds <- cbind(W = W, Z = Z, Y=Y)

  corrplot <- ggcorrplot(
    cor(ds),
    method = "square", 
    lab = TRUE, 
    lab_size = 3, 
    colors = c("#6D9EC1", "white", "#E46726"),
    title = "Correlation Matrix",
    ggtheme = theme_minimal()
  ) + 
    labs(fill = "Covariance")
}

