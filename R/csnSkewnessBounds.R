mu <- 0
sigma <- 1
gamma <- 10
nu <- 20
delta <- 1

# =========================================================================
# EMPIRICAL SKEWNESS COMPUTATION
# =========================================================================
n <- 10000
y <- csn::rcsn(k = n, mu = mu, sigma = sigma, gamma = gamma, nu = nu, delta = delta)

cat("=================================================================\n")
cat("Empirical Statistics\n")
cat("=================================================================\n")
cat(sprintf("Mean:     %.6f\n", mean(y)))
cat(sprintf("Std Dev:  %.6f\n", sd(y)))
cat(sprintf("Skewness: %.6f\n\n", moments::skewness(y)))

# =========================================================================
# THEORETICAL SKEWNESS COMPUTATION (from csn_1d_skewness.pdf)
# =========================================================================
S <- delta + gamma^2 * sigma
kappa <- nu / sqrt(S)
beta <- (gamma * sigma) / sqrt(S)
r0 <- dnorm(-kappa) / pnorm(-kappa)
kappa1 <- mu + beta * r0
kappa2 <- sigma + beta^2 * (kappa * r0 - r0^2)
kappa3 <- beta^3 * (r0 * (kappa^2 - 1) - 3 * kappa * r0^2 + 2 * r0^3)

cat("=================================================================\n")
cat("Theoretical Statistics\n")
cat("=================================================================\n")
cat(sprintf("Mean:     %.6f\n", kappa1))
cat(sprintf("Std Dev:  %.6f\n", sqrt(kappa2)))
cat(sprintf("Skewness: %.6f\n", kappa3 / (kappa2^(3/2))))


# =========================================================================
# THEORETICAL SKEWNESS CURVE OVER NU [-100, 100]
# =========================================================================
nu_grid <- -100:100

S <- delta + gamma^2 * sigma
beta <- (gamma * sigma) / sqrt(S)
kappa_vec <- nu_grid / sqrt(S)
r0_vec <- dnorm(kappa_vec) / pnorm(-kappa_vec)
kappa2_vec <- sigma + beta^2 * (kappa_vec * r0_vec - r0_vec^2)
kappa3_vec <- beta^3 * (r0_vec * (kappa_vec^2 - 1) - 3 * kappa_vec * r0_vec^2 + 2 * r0_vec^3)
skew_vec <- kappa3_vec / (kappa2_vec^(3/2))

plot(nu_grid, skew_vec,
     type = "l",
     xlab = "nu",
     ylab = "Theoretical skewness",
     main = "CSN theoretical skewness vs nu")
abline(h = 0, col = "gray70", lty = 2)

