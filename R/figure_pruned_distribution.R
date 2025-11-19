#' Illustration of pruned distribution
#'
#' This script replicates Figure 3 of the paper, comparing the pdf and cdf
#' of an univariate CSN distributed variable with its pruned approximation.
#'
#' @copyright 2022-2025 Gaygysyz Guljanov, Willi Mutschler, Mark Trede
#'
#' @note This is free software: you can redistribute it and/or modify
#' it under the terms of the GNU General Public License as published by
#' the Free Software Foundation, either version 3 of the License, or
#' (at your option) any later version.
#'
#' This is distributed in the hope that it will be useful,
#' but WITHOUT ANY WARRANTY; without even the implied warranty of
#' MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#' GNU General Public License for more details.
#'
#' This file is part of the replication files for the paper "Pruned Skewed
#' Kalman Filter and Smoother: With Application to DSGE models" by
#' Gaygysyz Guljanov, Willi Mutschler, and Mark Trede.

# housekeeping
setwd(this.path::here()) # set working directory
rm(list = ls()) # clear all variables and unload libraries

# load libraries
library(ggplot2)
library(csn)
library(patchwork)

# load pruned_csn_params function
source("csn_prune_params.R")

# define parameters
mu <- 0
sigma <- matrix(1, 1, 1)
gamma <- matrix(c(6, 0.1), 2, 1)
nu <- matrix(c(0, 0), 2, 1)
delta <- matrix(c(1, -0.1, -0.1, 1), 2, 2)

# compute unpruned parameters
cov_p <- rbind(cbind(sigma, sigma %*% t(gamma)), cbind(gamma %*% sigma, delta + gamma %*% sigma %*% t(gamma)))
corr_r <- cov2cor(cov_p)
cat("===================================\n")
cat("Correlation Matrix R Before Pruning\n")
cat("===================================\n")
print(round(corr_r, 4))

# compute pruned parameters
pruned_params <- pruned_csn_params(sigma, gamma, nu, delta, 0.1)
sigma_pruned <- pruned_params$sigma
gamma_pruned <- pruned_params$gamma
nu_pruned <- pruned_params$nu
delta_pruned <- pruned_params$delta

cat("Pruned parameters:\n")
cat("sigma_pruned =", as.numeric(sigma_pruned), "\n")
cat("gamma_pruned =", as.numeric(gamma_pruned), "\n")
cat("nu_pruned =", as.numeric(nu_pruned), "\n")
cat("delta_pruned =", as.numeric(delta_pruned), "\n\n")

# compute densities
x <- seq(-1, 4, length = 200)
densities <- data.frame(
  x = x,
  pdf_unpruned = dcsn(x, mu, sigma,        gamma,        as.vector(nu),        delta),
  pdf_pruned   = dcsn(x, mu, sigma_pruned, gamma_pruned, as.vector(nu_pruned), delta_pruned),
  cdf_unpruned = pcsn(x, mu, sigma,        gamma,        as.vector(nu),        delta),
  cdf_pruned   = pcsn(x, mu, sigma_pruned, gamma_pruned, as.vector(nu_pruned), delta_pruned)
)

# create plots
text_size <- 14 - 3
title_size <- 15 - 3
axis_text_size <- 12 - 3

p_pdfs <- ggplot(densities, aes(x = x)) +
  geom_line(aes(y = pdf_unpruned, linetype = "Unpruned"), linewidth = 1.2) +
  geom_line(aes(y = pdf_pruned, linetype = "Pruned"), linewidth = 1.2) +
  ggtitle("Probability density function") +
  ylab("") +
  scale_x_continuous(limits = c(-1, 4), expand = c(0, 0)) +
  scale_linetype_manual(
    name = "",
    values = c("Unpruned" = "solid", "Pruned" = "dashed")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    # Clean white background
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    # Add subtle panel border for better definition
    panel.border = element_rect(color = "grey85", fill = NA, linewidth = 0.5),
    # Grid styling for better readability
    panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    # Typography settings
    text = element_text(size = text_size),
    axis.title.y = element_text(size = text_size, margin = margin(r = 10)),
    axis.title.x = element_text(size = text_size, margin = margin(t = 10)),
    axis.text = element_text(size = axis_text_size, color = "grey30"),
    plot.title = element_text(size = title_size, hjust = 0.5, face = "plain", margin = margin(b = 10)),
    # Legend settings
    legend.title = element_blank(),
    legend.text = element_text(size = text_size),
    legend.key.width = unit(1.5, "cm"),
    legend.spacing.x = unit(0.5, "cm"),
    # Plot margins
    plot.margin = margin(10, 10, 10, 10)
  )
print(p_pdfs)

p_cdfs <- ggplot(densities, aes(x = x)) +
  geom_line(aes(y = cdf_unpruned, linetype = "Unpruned"), linewidth = 1.2) +
  geom_line(aes(y = cdf_pruned, linetype = "Pruned"), linewidth = 1.2) +
  ggtitle("Cumulative distribution function") +
  ylab("") +
  scale_x_continuous(limits = c(-1, 4), expand = c(0, 0)) +
  scale_linetype_manual(
    name = "",
    values = c("Unpruned" = "solid", "Pruned" = "dashed")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    # Clean white background
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    # Add subtle panel border for better definition
    panel.border = element_rect(color = "grey85", fill = NA, linewidth = 0.5),
    # Grid styling for better readability
    panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    # Typography settings
    text = element_text(size = text_size),
    axis.title.y = element_text(size = text_size, margin = margin(r = 10)),
    axis.title.x = element_text(size = text_size, margin = margin(t = 10)),
    axis.text = element_text(size = axis_text_size, color = "grey30"),
    plot.title = element_text(size = title_size, hjust = 0.5, face = "plain", margin = margin(b = 10)),
    # Legend settings
    legend.title = element_blank(),
    legend.text = element_text(size = text_size),
    legend.key.width = unit(1.5, "cm"),
    legend.spacing.x = unit(0.5, "cm"),
    # Plot margins
    plot.margin = margin(10, 10, 10, 10)
  )
print(p_cdfs)

fig_3 <- (p_pdfs + p_cdfs) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "center",
    legend.box.spacing = unit(0.1, "cm"),
    legend.margin = margin(t = 0, b = 5),
    plot.margin = margin(10, 10, 5, 10)
  )
print(fig_3)

# save plot
ggsave(
  filename = "figure_pruned_distribution.pdf",
  plot = fig_3,
  dpi = 300,
  units = "px",
  width = 3000,
  height = 1300
)

# =========================================================================
# COMPUTE KULLBACK-LEIBLER DIVERGENCE
# =========================================================================
# KL(P||Q) = ∫ p(x) log(p(x)/q(x)) dx, we'll compute this numerically using integration

# Determine integration range based on the distribution parameters
x_range_start <- min(-15, mu - 6 * sqrt(sigma[1, 1]))
x_range_end <- max(15, mu + 6 * sqrt(sigma[1, 1]))
x_kl <- seq(x_range_start, x_range_end, length = 20000)  # Use many points for better accuracy

# Compute densities
pdf_orig_kl <- dcsn(matrix(x_kl, length(x_kl), 1), mu, sigma, gamma, as.vector(nu), delta)
pdf_prun_kl <- dcsn(matrix(x_kl, length(x_kl), 1), mu, sigma_pruned, gamma_pruned, as.vector(nu_pruned), delta_pruned)

# Handle small densities carefully to avoid log(0)
eps <- 1e-20  # very small threshold
idx <- pdf_orig_kl > eps
x_kl_valid <- x_kl[idx]
pdf_orig_valid <- pdf_orig_kl[idx]
pdf_prun_valid <- pdf_prun_kl[idx]

# Compute integrand: p(x) * log(p(x)/q(x)) = p(x) * (log(p(x)) - log(q(x)))
log_ratio <- log(pdf_orig_valid) - log(pdf_prun_valid)
integrand <- pdf_orig_valid * log_ratio

# Numerical integration using trapezoidal rule
dx <- diff(x_kl_valid)
kl_divergence <- sum((integrand[-length(integrand)] + integrand[-1]) * dx / 2)

cat("==============================================================================\n")
cat("Computing Kullback-Leibler Divergence between original and pruned distribution\n")
cat("==============================================================================\n")
cat(sprintf("KL(P_original || P_pruned) = %.8f\n\n", kl_divergence))