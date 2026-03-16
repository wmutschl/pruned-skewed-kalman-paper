#' Prior and posterior density plots for shock parameters
#'
#' This script reads the kernel density estimates exported by Dynare
#' and creates faceted prior/posterior density plots for the standard
#' deviation and skewness parameters of the structural shocks.
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
library(dplyr)
library(ggplot2)
library(latex2exp)
library(patchwork)

# define architecture and MATLAB version
arch <- "maca64_m2max"
matlab_version <- "R2025b"

# load density data for both models
load_density <- function(model_name) {
  read.csv(
    paste0("../results/ireland2004/density_", model_name, "_", arch, "_", matlab_version, ".csv"),
    stringsAsFactors = FALSE
  ) %>%
    as_tibble() %>%
    mutate(model = model_name)
}

density_csn      <- load_density("csn")
density_gaussian <- load_density("gaussian")

# human-readable labels for parameters
param_labels <- c(
  "eta_a" = "Preference",
  "eta_e" = "Cost-Push",
  "eta_z" = "Productivity",
  "eta_r" = "Monetary Policy"
)

density_data <- bind_rows(density_csn, density_gaussian) %>%
  filter(type %in% c("stderr", "skew")) %>%
  mutate(
    parameter_label = factor(param_labels[parameter], levels = param_labels),
    distribution    = factor(distribution, levels = c("prior", "posterior"))
  )

# Okabe-Ito colorblind-safe palette
group_fills <- c(
  "Prior"                  = "#E69F00",
  "Posterior (CSN)"        = "#0072B2",
  "Posterior (Gaussian)"   = "#009E73"
)
group_outlines <- c(
  "Prior"                  = "#B47B00",
  "Posterior (CSN)"        = "#004D80",
  "Posterior (Gaussian)"   = "#006B5A"
)
group_linetypes <- c(
  "Prior"                  = "dashed",
  "Posterior (CSN)"        = "solid",
  "Posterior (Gaussian)"   = "dotdash"
)
group_alphas <- c(
  "Prior"                  = 0.25,
  "Posterior (CSN)"        = 0.50,
  "Posterior (Gaussian)"   = 0.35
)

# plot settings
text_size      <- 14
title_size     <- 16
axis_text_size <- 11
strip_size     <- 14

theme_paper <- theme_bw(base_size = 14) +
  theme(
    plot.background    = element_rect(fill = "white", color = NA),
    panel.background   = element_rect(fill = "white", color = NA),
    panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.6),
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    text               = element_text(size = text_size),
    axis.title.y       = element_text(size = text_size, margin = margin(r = 8)),
    axis.title.x       = element_text(size = text_size, margin = margin(t = 8)),
    axis.text          = element_text(size = axis_text_size, color = "black"),
    axis.ticks         = element_line(color = "black", linewidth = 0.4),
    strip.text         = element_text(size = strip_size, face = "bold"),
    strip.background   = element_rect(fill = "grey95", color = "grey70", linewidth = 0.4),
    legend.title       = element_blank(),
    legend.text        = element_text(size = text_size),
    legend.key.width   = unit(1.8, "cm"),
    legend.key.height  = unit(0.6, "cm"),
    legend.position    = "bottom",
    plot.margin        = margin(10, 12, 10, 10)
  )

make_density_plot <- function(data, x_label, plot_title, focus_posterior = FALSE,
                              padding_factor = 0.15) {
  all_levels <- c("Prior", "Posterior (CSN)", "Posterior (Gaussian)")
  data <- data %>%
    mutate(group = factor(group, levels = intersect(all_levels, unique(group))))
  lvls <- levels(data$group)

  p <- ggplot(data, aes(x = x, y = density, fill = group, color = group,
                        linetype = group, alpha = group)) +
    geom_area(linewidth = 0.7, position = "identity") +
    facet_wrap(~parameter_label, scales = "free", ncol = 4) +
    scale_fill_manual(values = group_fills[lvls]) +
    scale_color_manual(values = group_outlines[lvls]) +
    scale_linetype_manual(values = group_linetypes[lvls]) +
    scale_alpha_manual(values = group_alphas[lvls]) +
    labs(x = x_label, y = "Density", title = plot_title) +
    theme_paper +
    theme(plot.title = element_text(size = title_size, hjust = 0.5, face = "plain"))

  if (focus_posterior) {
    post_limits <- data %>%
      filter(group != "Prior") %>%
      group_by(parameter_label) %>%
      summarise(xmin = min(x), xmax = max(x), .groups = "drop") %>%
      mutate(pad = (xmax - xmin) * padding_factor,
             xmin = xmin - pad, xmax = xmax + pad)
    data <- data %>%
      left_join(post_limits, by = "parameter_label") %>%
      filter(x >= xmin, x <= xmax) %>%
      select(-xmin, -xmax, -pad)
    p <- ggplot(data, aes(x = x, y = density, fill = group, color = group,
                          linetype = group, alpha = group)) +
      geom_area(linewidth = 0.7, position = "identity") +
      facet_wrap(~parameter_label, scales = "free", ncol = 4) +
      scale_fill_manual(values = group_fills[lvls]) +
      scale_color_manual(values = group_outlines[lvls]) +
      scale_linetype_manual(values = group_linetypes[lvls]) +
      scale_alpha_manual(values = group_alphas[lvls]) +
      labs(x = x_label, y = "Density", title = plot_title) +
      theme_paper +
      theme(plot.title = element_text(size = title_size, hjust = 0.5, face = "plain"))
  }

  p
}

# --- Standard deviation parameters (Gaussian + CSN) ---
stderr_data <- density_data %>%
  filter(type == "stderr") %>%
  mutate(group = case_when(
    distribution == "prior" & model == "csn"          ~ "Prior",
    distribution == "posterior" & model == "csn"       ~ "Posterior (CSN)",
    distribution == "posterior" & model == "gaussian"  ~ "Posterior (Gaussian)"
  )) %>%
  filter(!is.na(group))

p_stderr <- make_density_plot(stderr_data, "Value", "Standard deviation of shocks",
                              focus_posterior = TRUE)

# --- Skewness parameters (CSN only) ---
skew_data <- density_data %>%
  filter(type == "skew", model == "csn") %>%
  mutate(group = if_else(distribution == "prior", "Prior", "Posterior (CSN)"))

p_skew <- make_density_plot(skew_data, "Value", "Skewness of shocks") +
  coord_cartesian(xlim = c(-1, 1))

# --- Combine into a single figure ---
fig_posterior <- p_stderr / p_skew +
  plot_layout(guides = "collect") &
  theme(
    legend.position      = "bottom",
    legend.box           = "horizontal",
    legend.justification = "center",
    legend.box.spacing   = unit(0.5, "cm"),
    legend.margin        = margin(t = 10, b = 5)
  )
print(fig_posterior)

# save figure
ggsave(
  filename = "figure_prior_posterior_plot.pdf",
  plot = fig_posterior,
  dpi = 300,
  units = "px",
  width = 3400,
  height = 2200
)
