get_density <- function(x, y, ...) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Please install MASS: install.packages('MASS')")
  }
  if (!is.numeric(x = x) | !is.numeric(x = y)) {
    stop("Must supply numeric values")
  }
  dens <- MASS::kde2d(x = x, y = y, ...)
  ix <- findInterval(x = x, vec = dens$x)
  iy <- findInterval(x = y, vec = dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

DensityScatter <- function(object, x, y, log_x = FALSE, log_y = FALSE, quantiles = NULL) {
  md <- object[[]]
  if (!(x %in% colnames(x = md))) {
    stop(x, " not found")
  }
  if (!(y %in% colnames(x = md))) {
    stop(y, " not found")
  }
  log10p <- function(x) {
    return(log10(x = x + 1))
  }
  null_fn <- function(x) {
    return(x)
  }
  logfnx <- ifelse(test = log_x, yes = log10p, no = null_fn)
  logfny <- ifelse(test = log_y, yes = log10p, no = null_fn)
  md$Density <- get_density(x = logfnx(md[[x]]), y = logfny(md[[y]]), n = 1000, h = c(1,
                                                                                      1))
  md <- md[order(md$Density), ]
  use_quantile <- FALSE
  if (!is.null(x = quantiles)) {
    use_quantile <- TRUE
    if (is.logical(x = quantiles)) {
      if (quantiles) {
        quantiles <- c(5, 10, 90, 95)
      } else {
        use_quantile <- FALSE
      }
    }
    if (!(all(quantiles >= 0) & all(quantiles <= 100))) {
      warning("Quantile values must be between 0 and 100", immediate. = TRUE)
      use_quantile <- FALSE
    }
    if (any(sapply(X = quantiles, FUN = function(x) x%%1 != 0))) {
      warning("Quantile values must be integers", immediate. = TRUE)
      use_quantile <- FALSE
    }
  }
  if (use_quantile) {
    quantiles <- paste0(quantiles, "%")
    x_quant <- quantile(x = md[[x]], probs = seq(0, 1, 0.01))
    y_quant <- quantile(x = md[[y]], probs = seq(0, 1, 0.01))
    xlines <- x_quant[quantiles]
    ylines <- y_quant[quantiles]
    xlines <- round(x = xlines, digits = 2)
    ylines <- round(x = ylines, digits = 2)
  }
  p <- ggplot(data = md, mapping = aes(x = .data[[x]], y = .data[[y]], color = .data[["Density"]])) +
    geom_point(size = 1) + scale_color_viridis_c(option = "B") + theme_bw()
  if (log_x) {
    p <- p + scale_x_log10()
  }
  if (log_y) {
    p <- p + scale_y_log10()
  }
  if (use_quantile) {
    p <- p + geom_vline(xintercept = unname(obj = xlines), color = "red") + 
      geom_hline(yintercept = unname(obj = ylines),color = "red") + 
      labs(title = "Quantiles", subtitle = paste0(x, ": ", paste0(names(x = xlines),
":", xlines, collapse = " "), "\n", y, ": ", paste0(names(x = ylines), ":", ylines,
collapse = " ")))
  }
  return(p)
}


