#' Estimated probability density functions of shocks
#'
#' This script replicates Figure 5 of the paper, comparing the
#' estimated probability density functions of shocks based on the
#' maximum likelihood estimates for the Gaussian and CSN models.
#'
#' @copyright 2026 Gaygysyz Guljanov, Willi Mutschler, Mark Trede
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
library(latex2exp)
library(patchwork)
utils::globalVariables(c("value")) # suppress R CMD check notes for ggplot2 non-standard evaluation

# define architecture and MATLAB version
arch <- "maca64_m4pro"
matlab_version <- "R2025b"

# load Gaussian estimates
mu_gauss <- as.matrix(read.csv(paste0("../results/ireland2004/ml/ireland2004_ml_1_gaussian_mu_e_", arch, "_", matlab_version, ".csv"), header = FALSE))
sigma_gauss <- as.matrix(read.csv(paste0("../results/ireland2004/ml/ireland2004_ml_1_gaussian_Sigma_e_", arch, "_", matlab_version, ".csv"), header = FALSE))
gamma_gauss <- as.matrix(read.csv(paste0("../results/ireland2004/ml/ireland2004_ml_1_gaussian_Gamma_e_", arch, "_", matlab_version, ".csv"), header = FALSE))
nu_gauss <- as.matrix(read.csv(paste0("../results/ireland2004/ml/ireland2004_ml_1_gaussian_nu_e_", arch, "_", matlab_version, ".csv"), header = FALSE))
delta_gauss <- as.matrix(read.csv(paste0("../results/ireland2004/ml/ireland2004_ml_1_gaussian_Delta_e_", arch, "_", matlab_version, ".csv"), header = FALSE))
smoothed_shocks_gauss <- read.csv(paste0("../results/ireland2004/ml/ireland2004_ml_1_gaussian_smoothed_shocks_", arch, "_", matlab_version, ".csv"))
eta_a_t_T_gauss <- data.frame(value = smoothed_shocks_gauss$eta_a) # smoothed normally distributed preference shock values
eta_e_t_T_gauss <- data.frame(value = smoothed_shocks_gauss$eta_e) # smoothed normally distributed cost-push shock values
eta_z_t_T_gauss <- data.frame(value = smoothed_shocks_gauss$eta_z) # smoothed normally distributed productivity shock values
eta_r_t_T_gauss <- data.frame(value = smoothed_shocks_gauss$eta_r) # smoothed normally distributed monetary policy shock values

# load CSN estimates
mu_csn <- as.matrix(read.csv(paste0("../results/ireland2004/ml/ireland2004_ml_3_csn_mu_e_", arch, "_", matlab_version, ".csv"), header = FALSE))
sigma_csn <- as.matrix(read.csv(paste0("../results/ireland2004/ml/ireland2004_ml_3_csn_Sigma_e_", arch, "_", matlab_version, ".csv"), header = FALSE))
gamma_csn <- as.matrix(read.csv(paste0("../results/ireland2004/ml/ireland2004_ml_3_csn_Gamma_e_", arch, "_", matlab_version, ".csv"), header = FALSE))
nu_csn <- as.matrix(read.csv(paste0("../results/ireland2004/ml/ireland2004_ml_3_csn_nu_e_", arch, "_", matlab_version, ".csv"), header = FALSE))
delta_csn <- as.matrix(read.csv(paste0("../results/ireland2004/ml/ireland2004_ml_3_csn_Delta_e_", arch, "_", matlab_version, ".csv"), header = FALSE))
smoothed_shocks_csn <- read.csv(paste0("../results/ireland2004/ml/ireland2004_ml_3_csn_smoothed_shocks_", arch, "_", matlab_version, ".csv"))
eta_a_t_T_csn <- data.frame(value = smoothed_shocks_csn$eta_a) # smoothed skew-normally distributed preference shock values
eta_e_t_T_csn <- data.frame(value = smoothed_shocks_csn$eta_e) # smoothed skew-normally distributed cost-push shock values
eta_z_t_T_csn <- data.frame(value = smoothed_shocks_csn$eta_z) # smoothed skew-normally distributed productivity shock values
eta_r_t_T_csn <- data.frame(value = smoothed_shocks_csn$eta_r) # smoothed skew-normally distributed monetary policy shock values

# define grid for plotting based on actual smoothed values with padding
# and compute optimal binwidth using Freedman-Diaconis rule
compute_grid <- function(gauss_vals, csn_vals, n_points = 200, padding_factor = 0.2,
                         binwidth_rule = c("fd", "scott", "sturges")) {
  binwidth_rule <- match.arg(binwidth_rule)
  all_vals <- c(gauss_vals, csn_vals)
  val_range <- range(all_vals)
  val_padding <- diff(val_range) * padding_factor
  x_grid <- seq(val_range[1] - val_padding, val_range[2] + val_padding, length = n_points)
  # Compute optimal binwidth based on selected rule
  binwidth <- switch(binwidth_rule,
    "fd" = 2 * IQR(all_vals) / (length(all_vals)^(1/3)),        # Freedman-Diaconis
    "scott" = 3.49 * sd(all_vals) / (length(all_vals)^(1/3)),   # Scott's rule
    "sturges" = diff(val_range) / ceiling(log2(length(all_vals)) + 1)  # Sturges (convert bins to width)
  )
  list(x = x_grid, binwidth = binwidth)
}

grid_eta_a <- compute_grid(eta_a_t_T_gauss$value, eta_a_t_T_csn$value, binwidth_rule = "scott")
grid_eta_e <- compute_grid(eta_e_t_T_gauss$value, eta_e_t_T_csn$value, binwidth_rule = "scott")
grid_eta_z <- compute_grid(eta_z_t_T_gauss$value, eta_z_t_T_csn$value, binwidth_rule = "scott")
grid_eta_r <- compute_grid(eta_r_t_T_gauss$value, eta_r_t_T_csn$value, binwidth_rule = "scott")

# compute densities
densities <- data.frame(
  x_eta_a = grid_eta_a$x,
  x_eta_e = grid_eta_e$x,
  x_eta_z = grid_eta_z$x,
  x_eta_r = grid_eta_r$x,
  f_eta_a_gauss = dnorm(grid_eta_a$x, 0, sqrt(sigma_gauss[1, 1])),
  f_eta_e_gauss = dnorm(grid_eta_e$x, 0, sqrt(sigma_gauss[2, 2])),
  f_eta_z_gauss = dnorm(grid_eta_z$x, 0, sqrt(sigma_gauss[3, 3])),
  f_eta_r_gauss = dnorm(grid_eta_r$x, 0, sqrt(sigma_gauss[4, 4])),
  f_eta_a_csn = dcsn(grid_eta_a$x, mu_csn[1, 1], sigma_csn[1, 1], gamma_csn[1, 1], nu_csn[1, 1], delta_csn[1, 1]),
  f_eta_e_csn = dcsn(grid_eta_e$x, mu_csn[2, 1], sigma_csn[2, 2], gamma_csn[2, 2], nu_csn[2, 1], delta_csn[2, 2]),
  f_eta_z_csn = dcsn(grid_eta_z$x, mu_csn[3, 1], sigma_csn[3, 3], gamma_csn[3, 3], nu_csn[3, 1], delta_csn[3, 3]),
  f_eta_r_csn = dcsn(grid_eta_r$x, mu_csn[4, 1], sigma_csn[4, 4], gamma_csn[4, 4], nu_csn[4, 1], delta_csn[4, 4])
)

# print optimal binwidths for reference
cat("Optimal binwidths (Scott's rule):\n")
cat("  eta_a:", round(grid_eta_a$binwidth, 4), "\n")
cat("  eta_e:", round(grid_eta_e$binwidth, 4), "\n")
cat("  eta_z:", round(grid_eta_z$binwidth, 4), "\n")
cat("  eta_r:", round(grid_eta_r$binwidth, 4), "\n")

# plot settings
text_size <- 14
title_size <- 15 - 1
axis_text_size <- 12 - 3

# preference shock plot
create_shock_plot <- function(csn_data, gauss_data, grid_info, densities, x_col, f_csn_col, f_gauss_col,
                              x_label, plot_title, show_legend = FALSE, show_y_label = TRUE) {
  p <- ggplot(data = densities, aes(x = .data[[x_col]])) +
    geom_histogram(data = csn_data, aes(x = .data[["value"]], y = after_stat(density)), inherit.aes = FALSE,
                   binwidth = grid_info$binwidth, alpha = 0.4, color = "white", fill = "steelblue",
                   linewidth = 0.3) +
    geom_histogram(data = gauss_data, aes(x = .data[["value"]], y = after_stat(density)), inherit.aes = FALSE,
                   binwidth = grid_info$binwidth, alpha = 0.4, color = "white", fill = "darkorange",
                   linewidth = 0.3) +
    geom_line(aes(y = .data[[f_csn_col]], color = "Skew-Normal"), linewidth = 1.5, linetype = "solid") +
    geom_line(aes(y = .data[[f_gauss_col]], color = "Gaussian"), linewidth = 1.5, linetype = "dashed") +
    scale_color_manual(values = c("Skew-Normal" = "steelblue", "Gaussian" = "darkorange")) +
    labs(
      x = x_label,
      y = if(show_y_label) "Density" else "",
      title = plot_title
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
      plot.title = element_text(size = title_size, hjust = 0.5, face = "plain",
                                margin = margin(b = 10)),
      # Legend settings
      legend.title = element_blank(),
      legend.text = element_text(size = text_size),
      legend.key.width = unit(1.5, "cm"),
      legend.spacing.x = unit(0.5, "cm"),
      # Plot margins
      plot.margin = margin(10, 10, 10, 10)
    )
  return(p)
}

# preference shock plot
p_pref <- create_shock_plot(
  csn_data = eta_a_t_T_csn,
  gauss_data = eta_a_t_T_gauss,
  grid_info = grid_eta_a,
  densities = densities,
  x_col = "x_eta_a",
  f_csn_col = "f_eta_a_csn",
  f_gauss_col = "f_eta_a_gauss",
  x_label = TeX(r"($\eta_{a}$)"),
  plot_title = "Preference shock",
  show_legend = TRUE
)
print(p_pref)

# cost-push shock plot
p_cost <- create_shock_plot(
  csn_data = eta_e_t_T_csn,
  gauss_data = eta_e_t_T_gauss,
  grid_info = grid_eta_e,
  densities = densities,
  x_col = "x_eta_e",
  f_csn_col = "f_eta_e_csn",
  f_gauss_col = "f_eta_e_gauss",
  x_label = TeX(r"($\eta_{e}$)"),
  plot_title = "Cost-push shock",
  show_legend = TRUE,
  show_y_label = FALSE
)
print(p_cost)

# productivity shock plot
p_prod <- create_shock_plot(
  csn_data = eta_z_t_T_csn,
  gauss_data = eta_z_t_T_gauss,
  grid_info = grid_eta_z,
  densities = densities,
  x_col = "x_eta_z",
  f_csn_col = "f_eta_z_csn",
  f_gauss_col = "f_eta_z_gauss",
  x_label = TeX(r"($\eta_{z}$)"),
  plot_title = "Productivity shock",
  show_legend = TRUE
)
print(p_prod)

# monetary policy shock plot
p_monp <- create_shock_plot(
  csn_data = eta_r_t_T_csn,
  gauss_data = eta_r_t_T_gauss,
  grid_info = grid_eta_r,
  densities = densities,
  x_col = "x_eta_r",
  f_csn_col = "f_eta_r_csn",
  f_gauss_col = "f_eta_r_gauss",
  x_label = TeX(r"($\eta_{r}$)"),
  plot_title = "Monetary policy shock",
  show_legend = TRUE,
  show_y_label = FALSE
)
print(p_monp)

# create figure
fig_4 <- (p_pref + p_cost) / (p_prod + p_monp) +
  plot_layout(guides = "collect", heights = c(1, 1)) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "center",
    legend.box.spacing = unit(0.5, "cm"),
    legend.margin = margin(t = 10, b = 5),
    plot.margin = margin(10, 10, 10, 10)
  )
print(fig_4)

# save figure
ggsave(
  filename = "figure_estimated_distributions.pdf",
  plot = fig_4,
  dpi = 300,
  units = "px",
  width = 3200,
  height = 2800
)
