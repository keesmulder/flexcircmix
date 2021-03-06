
#' Circular ggplot scales
#'
#' Extensions of \code{scale_x_continuous} and \code{scale_y_continuous} that
#' are useful for circular axes.
#'
#' Behind the screens, radians will be used, but the plot will display
#' transformed values.
#'
#' The continuous transformations are \code{"radians"}, \code{"degrees"} and
#' \code{"hours"}. For these, \code{nticks} and \code{digits} can be set.
#'
#' The categorical transformations are \code{"cardinal"}, which prints character
#' labels of the directions (left, right, up, down), and \code{"compass"}, which
#' prints compass directions (North, South, East, West). Both of these assume
#' the circular data starts on the right and moves counterclockwise. If not,
#' appropriate transformations must first be taken, or the labels will be
#' nonsensical.
#'
#' Two additional labels are \code{"texpi"} and \code{"texnegpi"}, which print
#' nice labels for TeX to interpret starting at 0 and \code{-pi} respectively.
#' This is useful if TikZ is used.
#'
#' @param units The units to display on the axis. See `Details` for the options.
#' @param nticks The number of ticks to display. Only relevant for continuous
#'   scales.
#' @param digits Numeric; The number of digits for continuous scales.
#' @param limits Numeric vector; The limits of the plot. Must two numbers that
#'   are  \code{2*pi} apart.
#' @param scale_function The scale function from ggplot to use. Either
#'   \code{scale_x_continuous} or \code{scale_y_continuous}. For convenience
#'   this option can be circumvented by directly using \code{scale_x_circular}
#'   or \code{scale_y_circular}.
#' @param ... Additional arguments to be passed to  \code{scale_x_continuous} or
#'   \code{scale_y_continuous}.
#'
#' @return An object of type \code{ScaleContinuousPosition} that can be added to
#'   any existing \code{ggplot}.
#' @export
#'
#' @examples
#' p <- plot_batmixfit(params = cbind(mu = c(-pi/2, 0, pi/2, pi), kp = 4,
#'                               lam = c(-.9, .2, .8, 0), alph = .25))
#'
#' p + scale_x_circular(units = "compass")
#' p + scale_x_circular(units = "compass", limits = c(-1, 2*pi - 1))
#'
#' p + scale_x_circular(units = "cardinal")
#'
#' p + scale_x_circular(units = "texpi")
#' p + scale_x_circular(units = "texnegpi")
#'
#' p + scale_x_circular(units = "hours", nticks = 24)
#'
scale_circular <- function(units = "degrees", nticks = 4,
                           digits = 0, limits = c(0, 2 * pi),
                           scale_function = ggplot2::scale_x_continuous, ...) {

  if (round(abs(limits[1] - limits[2]) - 2*pi, 3) != 0) {
    stop("Limits must be a range of length 2*pi.")
  }

  units_categ <- units %in% c("texpi", "texnegpi", "cardinal", "compass")

  if (units_categ & nticks != 4) {
    warning(paste0("For `units = ", units, "`, setting `nticks` to 4."))
    nticks <- 4
  }

  if (units_categ) {
    if (units == "texpi") {
      brks   <- seq(0, 2 * pi, length.out = 5)
      limits <- c(0, 2 * pi)
    } else if (units == "texnegpi") {
      brks   <- seq(-pi, pi, length.out = 5)
      limits <- c(-pi, pi)
    } else {
      # These are the options for cardinal or compass directions.
      brks   <- seq(-pi, 3*pi, length.out = 9)
    }
  } else {
    brks <- seq(limits[1], limits[2], length.out = nticks + 1)
  }

  conv_brks <- switch(
    units,
    radians  = round(brks, digits),
    degrees  = round(brks * 180 / pi, digits),
    hours    = round(brks * 12 / pi, digits),
    texpi    = c("$0$", "$\\frac{\\pi}{2}$", "$\\pi$",
                 "$\\frac{3\\pi}{2}$", "$2\\pi$"),
    texnegpi = c("$-\\pi$", "$-\\frac{\\pi}{2}$", "$0$",
                 "$\\frac{\\pi}{2}$", "$\\pi$"),
    cardinal = c("Right", "Up", "Left", "Down")[c(3:4, 1:4, 1:3)],
    compass  = c("East", "North", "West", "South")[c(3:4, 1:4, 1:3)])

  scale_function(breaks = brks,
                 labels = conv_brks,
                 limits = limits,
                 ...)
}

#' @export
#' @describeIn scale_circular X axis has circular scale.
scale_x_circular <- function(...) {
  scale_circular(scale_function = ggplot2::scale_x_continuous, ...)
}

#' @export
#' @describeIn scale_circular Y axis has circular scale.
scale_y_circular <- function(...) {
  scale_circular(scale_function = ggplot2::scale_y_continuous, ...)
}




#' Plot a mixture of Batschelet type distributions.
#'
#' @param x Optional data vector of angles in radians to plot a histogram of.
#' @param params A matrix of parameters.
#' @param dbat_fun A function; the Batschelet function to use. Defaults to the
#'   Inverse Batschelet Function.
#' @param bins Integer; The number of bins to use in the optional histogram.
#' @param hist_alpha Numeric; The alpha value of the histogram of the
#'   data.
#' @param res Integer; The number of points at which to evaluate the density
#'   function.
#'
#' @return A ggplot.
#' @export
#'
#' @examples
#' plot_batmixfit(params = cbind(mu = c(-pi/2, 0, pi/2, pi), kp = 4,
#'                               lam = c(-.9, .2, .8, 0), alph = .25))
#'
plot_batmixfit <- function(x, params, dbat_fun = dpowbat, bins = 100, res = 400,
                           hist_alpha = .3) {


  # Initialize plot.
  if (missing(x)) {
    p <- ggplot2::ggplot(data.frame(x = c(-pi, pi)))
  } else {
    p <- ggplot2::ggplot(data.frame(x = force_neg_pi_pi(x)))
  }

  p <- p +
    ggplot2::xlim(-pi, pi) +
    ggplot2::ylim(c(0, NA)) +
    ggplot2::coord_cartesian(expand = TRUE) +
    ggplot2::theme_bw()


  # If we don't have data given as well, return only the pdfs. Otherwise, add a
  # histogram.
  if (!missing(x)) {
    p <- p +
      ggplot2::geom_histogram(
        mapping = ggplot2::aes_string(x = "x", y = "..density.."),
        fill = grDevices::rgb(.65, .65, .85, hist_alpha), col = "white",
        boundary = -pi, binwidth = 2*pi / bins)
  }


  # Add the mixture density.
  p <- p +
    ggplot2::stat_function(fun = dbatmix_pmat,
                           args = list(pmat = params, dbat_fun = dbat_fun),
                           n = res,
                           colour = grDevices::rgb(.35, .35, .35, .8),
                           lwd = 1.2)


  # Add the separate densities.
  n_comp <- nrow(params)

  # Force evaluation of the functions parameters, because otherwise lazy
  # evaluation will only cause us to plot the last function.
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
#' @param ... Additional arguments to be passed to \code{plot_batmixfit}.
#'
#' @return A ggplot.
#' @export
#'
#' @examples
#'
#' \dontrun{
#'   movMF_model <- movMF(cbind(runif(100, -1, 1), runif(100, -1, 1), 2), k = 4)
#'   plot_movMF_as_batmix(movMF_model)
#' }
#'
plot_movMF_as_batmix <- function(m, ...) {
  pmat <- invbatmix_pmat_from_movMF(m)
  plot_batmixfit(params = pmat, ...)
}




#' Plot a Batschelet-type mixture model
#'
#' @param x A \code{batmixmod} object.
#' @param ... Additional arguments to be passed to \code{plot_batmixfit}.
#'
#' @return A \code{ggplot}.
#' @export
#'
#' @examples
#' x <- rinvbatmix(50)
#' plot(fitbatmix(x, method = "EM"))
#'
plot.batmixmod <- function(x, ...) {
  plot_batmixfit(x$x, params = x$estimates, ...)
}



#' Plot a sample of Batschelet mixture parameter sets.
#'
#' @param x An optional dataset of angles to be plotted as a histogram.
#' @param param A matrix of parameter sets.
#' @param dbat_fun A function; The pdf of the chosen Batschelet distribution.
#' @param bins The number of bins to draw in the histogram.
#' @param res Number of points at which to evaluate the functions.
#' @param orderColor Logical; If \code{TRUE}, plotted pdfs get darker and more
#'   red as they were sampled later.
#' @param plot_n Integer; the number of parameter rows to sample from the matrix
#'   \code{param}. This is intended for MCMC for example, where we can take a
#'   subsample of the parameter matrix to plot for speed.
#' @param hist_alpha Numeric; The alpha value of the histogram of the
#'   data.
#' @param dens_darkness Numeric; Higher numbers result in less transparent
#'   densities plotted.
#'
#' @return A ggplot.
#' @export
#'
#' @examples
#' x <- rinvbatmix(50)
#' mod <- fitbatmix(x, method = "bayes",  Q = 10)
#'
#' plot_batmix_sample(x, mod$mcmc_sample, dens_darkness = 5)
#'
plot_batmix_sample <- function(x, param, dbat_fun = dinvbat,
                               plot_n = nrow(param),
                               hist_alpha = .3, dens_darkness = 20,
                               bins = 100, res = 400, orderColor = FALSE) {

  # Change to matrix if needed.
  if (is.vector(param)) param <- t(param)

  # Initialize plot.
  if (missing(x)) {
    p <- ggplot2::ggplot(data.frame(x = c(-pi, pi)))
  } else {
    p <- ggplot2::ggplot(data.frame(x))
  }

  p <- p +
    ggplot2::xlim(-pi, pi) +
    ggplot2::ylim(c(0, NA)) +
    ggplot2::coord_cartesian(expand = TRUE) +
    ggplot2::theme_bw()


  # If we don't have data given as well, return only the pdfs. Otherwise, add a
  # histogram.
  if (!missing(x)) {
    p <- p +
      ggplot2::geom_histogram(mapping = ggplot2::aes_string(x = "x",
                                                            y = "..density.."),
                              fill = grDevices::rgb(.65, .65, .85, hist_alpha),
                              col = "white",
                              boundary = -pi, binwidth = 2*pi / bins)

    # If there are no parameters given, just return the histogram.
    if (identical(param, NA)) return(p)
  }


  # Remove some rows if we don't plot every row of param.
  param <- param[round(seq(1, nrow(param),
                           length.out = plot_n)), , drop = FALSE]

  mu_mat   <- param[, grep("mu_[0-9]",   colnames(param)), drop = FALSE]
  kp_mat   <- param[, grep("kp_[0-9]",   colnames(param)), drop = FALSE]
  lam_mat  <- param[, grep("lam_[0-9]",  colnames(param)), drop = FALSE]
  alph_mat <- param[, grep("alph_[0-9]", colnames(param)), drop = FALSE]

  n_comp <- ncol(mu_mat)

  if (orderColor) ordseq <- seq(0, 1, 0.6/plot_n)


  for (i in 1:plot_n) {
    suppressWarnings(
      p <- p + ggplot2::stat_function(
        fun = dbatmix,
        args = list(mus      = mu_mat[i, ],
                    kps      = kp_mat[i, ],
                    lams     = lam_mat[i, ],
                    alphs    = alph_mat[i, ],
                    dbat_fun = dbat_fun),
        col = grDevices::rgb(ifelse(orderColor, ordseq[i], 0.2),
                             0.2, 0.2, min(1, dens_darkness/plot_n)),
        n = res)
    )
  }
  p
}


