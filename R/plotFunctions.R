plot_invbatmixfit <- function(x, params, bins = 100, res = 400) {

  # Initialize plot.
  p <- ggplot2::ggplot(data.frame(x)) +
    ggplot2::xlim(-pi, pi) +
    ggplot2::ylim(c(0, NA)) +
    ggplot2::coord_cartesian(expand = TRUE) +
    ggplot2::theme_bw()


  # If we don't have data given as well, return only the pdfs. Otherwise, add a histogram.
  if (!missing(x)) {
    p <- p +
      ggplot2::geom_histogram(mapping = ggplot2::aes_string(x = "x", y = "..density.."),
                              fill = rgb(.65, .65, .85, .3), col = "white",
                              boundary = -pi, binwidth = 2*pi / bins)
  }

  # Add the mixture density.
  p <- p +
    ggplot2::stat_function(fun = dinvbatmix_pmat,
                           args = list(pmat = params), n = res,
                           colour = rgb(.35, .35, .35, .8),
                           lwd = 1.2)


  # Add the separate densities.
  n_comp <- nrow(params)

  # Force evaluation of the functions parameters, because otherwise lazy evaluation will only cause
  # us to plot the last function.
  funlist <- lapply(1:n_comp, function(compi) {
    function(x) params[compi, 'alph'] * dinvbat(x,
                                                mu = params[compi, 'mu'],
                                                kp = params[compi, 'kp'],
                                                lam = params[compi, 'lam'])})

  if (n_comp < 10) {
    palette <- RColorBrewer::brewer.pal(n_comp, "Set1")
  } else {
    palette <- RColorBrewer::brewer.pal(n_comp, "RdYlBu")
  }

  # Next, add each separate pdf.
  for (compi in 1:n_comp) {
    p <- p + ggplot2::stat_function(fun = funlist[[compi]],
                                    col = palette[compi],
                                    lwd = .8, n = res)
  }
  p
}


plot_movMF_as_invbatmix <- function(x, m, bins = 100, res = 400) {
  pmat <- invbatmix_pmat_from_movMF(m)
  plot_invbatmixfit(x, pmat, bins = bins, res = res)
}



