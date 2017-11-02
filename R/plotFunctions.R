#' Plot a mixture of Batschelet type distributions.
#'
#' @param x Optional data vector of angles in radians to plot a histogram of.
#' @param params A matrix of parameters.
#' @param dbat_fun A function; the Batschelet function to use. Defaults to the Inverse Batschelet Function.
#' @param bins Integer; The number of bins to use in the optional histogram.
#' @param res Integer; The number of points at which to evaluate the density function.
#'
#' @return A ggplot.
#' @export
#'
#' @examples
#' plot_batmixfit(params = cbind(mu = c(-pi/2, 0, pi/2, pi), kp = 4, lam = c(-.9, .2, .8, 0), alph = .25))
#'
plot_batmixfit <- function(x, params, dbat_fun = dinvbat, bins = 100, res = 400) {


  # Initialize plot.
  if (missing(x)) {
    p <- ggplot2::ggplot(data.frame(x = c(-pi, pi)))
  } else {
    p <- ggplot2::ggplot(data.frame(x))
  }

  p <-p +
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
    ggplot2::stat_function(fun = dbatmix_pmat,
                           args = list(pmat = params, dbat_fun = dbat_fun),
                           n = res,
                           colour = rgb(.35, .35, .35, .8),
                           lwd = 1.2)


  # Add the separate densities.
  n_comp <- nrow(params)

  # Force evaluation of the functions parameters, because otherwise lazy evaluation will only cause
  # us to plot the last function.
  funlist <- lapply(1:n_comp, function(compi) {
    function(x) params[compi, 'alph'] * dbat_fun(x,
                                                mu = params[compi, 'mu'],
                                                kp = params[compi, 'kp'],
                                                lam = params[compi, 'lam'])})

  if (n_comp < 10) {
    palette <- RColorBrewer::brewer.pal(max(n_comp, 3), "Set1")
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


#' plot_movMF_as_batmix
#'
#' Plot a result from movMF using the Batschelet mixture plotting functions.
#'
#' @param m A movMF results object.
#'
#' @return A ggplot.
#' @export
#'
#' @examples
#' movMF_model <- movMF(cbind(runif(100, -1, 1), runif(100, -1, 1), 2), k = 4)
#' plot_movMF_as_batmix(movMF_model)
#'
plot_movMF_as_batmix <- function(m, ...) {
  pmat <- invbatmix_pmat_from_movMF(m)
  plot_batmixfit(params = pmat, ...)
}



