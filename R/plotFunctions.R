plot_invbatmixfit <- function(x, params) {
  ggplot2::ggplot(data.frame(x)) +
    ggplot2::geom_histogram(ggplot2::aes_string(x = "x", y = "..density.."), binwidth = 4*pi/400) +
    ggplot2::stat_function(fun = dinvbatmix,
                           args = list(mu = params[, 'mu'], kp = params[, 'kp'],
                                       lam = params[, 'lam'], alph = params[, 'alph']),
                           lwd = 1.2, col = "skyblue") +
    ggplot2::xlim(-pi - .1, pi + .1) + ggplot2::theme_bw()
}