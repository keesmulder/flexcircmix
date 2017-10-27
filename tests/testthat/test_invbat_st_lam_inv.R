library(flexcircmix)

context("T and S of lambda and their inverses")


curve(s_lam(x, .8), -pi, pi, asp = 1)
curve(s_lam(x, -.8), -pi, pi, asp = 1, add = TRUE)
curve(s_lam(x, 1), -pi, pi, asp = 1, add = TRUE)
curve(s_lam(x, -1), -pi, pi, asp = 1, add = TRUE)
curve(s_lam(x, 0), -pi, pi, asp = 1, add = TRUE)
# curve(s_lam_approx(x, .8), -pi, pi, add = TRUE, col = "skyblue")



curve(t_lam(x, .8), -pi, pi, asp = 1)
curve(t_lam(x, -.8), -pi, pi, asp = 1, add = TRUE)
curve(t_lam(x, 1), -pi, pi, asp = 1, add = TRUE)
curve(t_lam(x, -1), -pi, pi, asp = 1, add = TRUE)
curve(t_lam(x, 0), -pi, pi, asp = 1, add = TRUE)




# curve(tpow_lam(x, -1.02), -pi, pi, asp = 1)
# curve(tpow_lam_d(x, -1.02), -pi, pi, asp = 1, add = TRUE)
#
# curve(tpow_lam_d(x, -.01), -pi, pi, asp = 1)
# curve(tpow_lam(x, -.01), -pi, pi, asp = 1, add = TRUE)


th <- force_neg_pi_pi(circular:::RvonmisesRad(100000, 0, 1))

th2 <- tpow_lam(th, lam = .2)

plot(density(th))
plot(density(th2))
curve(dpowbat(x, kp = 1, lam = .2), -pi, pi, add = TRUE, col = "goldenrod")
curve(dpowbat(x, kp = 1, lam = 1), -pi, pi, add = TRUE, col = "goldenrod")
curve(dpowbat(x, kp = 1, lam = 2), -pi, pi, add = TRUE, col = "goldenrod")
curve(dpowbat(x, kp = 1, lam = 0), -pi, pi, add = TRUE, col = "goldenrod")

curve(dpowbat(x, kp = 1, lam = 1), -pi, pi, add = TRUE, col = "goldenrod")



curve(exp(cos(x - sin(x))), -pi, pi)
curve(exp(cos(x + sin(x))), -pi, pi)

curve(x + 1 * sin(x), -pi, pi, asp = 1)
curve(x - sin(x), -pi, pi, add = TRUE, col = "skyblue")
curve(x + 0, -pi, pi, add = TRUE, col = "tomato")


# Some good options for approximations
curve(sign(x) * pi*(abs(x)/pi)^(1/3), -pi, pi, add = TRUE, col = "goldenrod")
curve(sign(x) * pi*(abs(x)/pi)^(1/2.3), -pi, pi, add = TRUE, col = "goldenrod")


# Original functions
curve(t_lam(x, 1), -pi, pi, asp = 1, lty = "dashed")
curve(t_lam(x, -1), -pi, pi, asp = 1, add = TRUE, lty = "dashed")
curve(t_lam(x, .5), -pi, pi, asp = 1, add = TRUE, lty = "dotted")
curve(t_lam(x, -.5), -pi, pi, asp = 1, add = TRUE, lty = "dotted")
curve(t_lam(x, 0), -pi, pi, asp = 1, add = TRUE)

# Some good options for approximations
curve(sign(x) * pi*(abs(x)/pi)^(1/3), -pi, pi, add = TRUE, col = "goldenrod")
curve(sign(x) * pi*(abs(x)/pi)^(3), -pi, pi, add = TRUE, col = "goldenrod")
curve(sign(x) * pi*(abs(x)/pi)^(1/2), -pi, pi, add = TRUE, col = "goldenrod")
curve(sign(x) * pi*(abs(x)/pi)^(2), -pi, pi, add = TRUE, col = "goldenrod")
curve(sign(x) * pi*(abs(x)/pi)^(1/sqrt(2)), -pi, pi, add = TRUE, col = "goldenrod")
curve(sign(x) * pi*(abs(x)/pi)^(sqrt(2)), -pi, pi, add = TRUE, col = "goldenrod")



