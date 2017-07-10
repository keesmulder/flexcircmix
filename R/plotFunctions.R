plot_data_fit <- function(x, param) {
  ggplot2::ggplot(data.frame(x)) +
    ggplot2::geom_histogram(ggplot2::aes_string(x = "x", y = "..density.."), binwidth = 4*pi/400) +
    ggplot2::stat_function(fun = dinvbat, args = list(mu = param[1], kp = param[2], lam = param[3])) +
    ggplot2::xlim(-pi, pi) + ggplot2::theme_bw()
}