#' Prior and posterior density plots for shock and deep parameters
#'
#' This script reads the kernel density estimates exported by Dynare
#' and creates faceted prior/posterior density plots for the standard
#' deviation, skewness, and deep structural parameters.
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
arch <- "maca64_m4pro"
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

# human-readable labels for shock parameters
shock_labels <- c(
  "eta_a" = "Preference",
  "eta_e" = "Cost-Push",
  "eta_z" = "Productivity",
  "eta_r" = "Monetary Policy"
)

# human-readable labels for deep parameters (plotmath strings for label_parsed)
deep_labels <- c(
  "OMEGA"  = "omega",
  "RHO_PI" = "rho[pi]",
  "RHO_G"  = "rho[g]",
  "RHO_X"  = "rho[x]",
  "RHO_A"  = "rho[a]",
  "RHO_E"  = "rho[e]"
)

density_data <- bind_rows(density_csn, density_gaussian) %>%
  filter(type %in% c("stderr", "skew", "deep")) %>%
  mutate(
    parameter_label = case_when(
      type %in% c("stderr", "skew") ~ shock_labels[parameter],
      type == "deep"                 ~ parameter
    ),
    distribution = factor(distribution, levels = c("prior", "posterior"))
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
title_size     <- 15
axis_text_size <- 12
strip_size     <- 13

theme_paper <- theme_minimal(base_size = 14) +
  theme(
    plot.background    = element_rect(fill = "white", color = NA),
    panel.background   = element_rect(fill = "white", color = NA),
    panel.border       = element_rect(color = "grey85", fill = NA, linewidth = 0.5),
    panel.grid.major   = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.minor   = element_blank(),
    text               = element_text(size = text_size),
    axis.title.y       = element_text(size = text_size, margin = margin(r = 10)),
    axis.title.x       = element_text(size = text_size, margin = margin(t = 10)),
    axis.text          = element_text(size = axis_text_size, color = "grey30"),
    strip.text         = element_text(size = strip_size),
    strip.background   = element_blank(),
    legend.title       = element_blank(),
    legend.text        = element_text(size = text_size),
    legend.key.width   = unit(1.5, "cm"),
    legend.spacing.x   = unit(0.5, "cm"),
    legend.position    = "bottom",
    plot.margin        = margin(10, 10, 10, 10)
  )

make_density_plot <- function(data, x_label, plot_title, focus_posterior = FALSE,
                              padding_factor = 0.15, ncol = 4, facet_labeller = label_value) {
  all_levels <- c("Prior", "Posterior (CSN)", "Posterior (Gaussian)")
  data <- data %>%
    mutate(group = factor(group, levels = intersect(all_levels, unique(group))))
  lvls <- levels(data$group)

  build_plot <- function(d) {
    ggplot(d, aes(x = x, y = density, fill = group, color = group,
                  linetype = group, alpha = group)) +
      geom_area(linewidth = 0.7, position = "identity") +
      facet_wrap(~parameter_label, scales = "free", ncol = ncol,
                 labeller = facet_labeller) +
      scale_fill_manual(values = group_fills[lvls]) +
      scale_color_manual(values = group_outlines[lvls]) +
      scale_linetype_manual(values = group_linetypes[lvls]) +
      scale_alpha_manual(values = group_alphas[lvls]) +
      labs(x = x_label, y = "Density", title = plot_title) +
      theme_paper +
      theme(plot.title = element_text(size = title_size, hjust = 0.5, face = "plain",
                                      margin = margin(b = 10)))
  }

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
  }

  build_plot(data)
}

# --- Standard deviation parameters (Gaussian + CSN) ---
stderr_data <- density_data %>%
  filter(type == "stderr") %>%
  mutate(
    parameter_label = factor(parameter_label, levels = shock_labels),
    group = case_when(
      distribution == "prior" & model == "csn"          ~ "Prior",
      distribution == "posterior" & model == "csn"       ~ "Posterior (CSN)",
      distribution == "posterior" & model == "gaussian"  ~ "Posterior (Gaussian)"
    )
  ) %>%
  filter(!is.na(group))

p_stderr <- make_density_plot(stderr_data, "Value", "Standard deviation of shocks",
                              focus_posterior = TRUE)

# --- Skewness parameters (CSN only) ---
skew_data <- density_data %>%
  filter(type == "skew", model == "csn") %>%
  mutate(
    parameter_label = factor(parameter_label, levels = shock_labels),
    group = if_else(distribution == "prior", "Prior", "Posterior (CSN)")
  )

p_skew <- make_density_plot(skew_data, "Value", "Skewness of shocks") +
  coord_cartesian(xlim = c(-1, 1))

# --- Combine shock parameters into a single figure ---
fig_shocks <- p_stderr / p_skew +
  plot_layout(guides = "collect") &
  theme(
    legend.position      = "bottom",
    legend.box           = "horizontal",
    legend.justification = "center",
    legend.box.spacing   = unit(0.5, "cm"),
    legend.margin        = margin(t = 10, b = 5),
    plot.margin          = margin(10, 10, 10, 10)
  )
print(fig_shocks)

ggsave(
  filename = "figure_prior_posterior_shocks.pdf",
  plot = fig_shocks,
  dpi = 300,
  units = "px",
  width = 3400,
  height = 2200
)

# --- Deep parameters (Gaussian + CSN) ---
deep_param_labeller <- as_labeller(deep_labels, default = label_parsed)

deep_data <- density_data %>%
  filter(type == "deep") %>%
  mutate(
    parameter_label = factor(parameter_label, levels = names(deep_labels)),
    group = case_when(
      distribution == "prior" & model == "csn"          ~ "Prior",
      distribution == "posterior" & model == "csn"       ~ "Posterior (CSN)",
      distribution == "posterior" & model == "gaussian"  ~ "Posterior (Gaussian)"
    )
  ) %>%
  filter(!is.na(group))

p_deep <- make_density_plot(deep_data, "Value", "Deep parameters",
                            focus_posterior = TRUE, ncol = 3,
                            facet_labeller = deep_param_labeller)

fig_deep <- p_deep +
  theme(
    legend.position      = "bottom",
    legend.box           = "horizontal",
    legend.justification = "center",
    legend.box.spacing   = unit(0.5, "cm"),
    legend.margin        = margin(t = 10, b = 5),
    plot.margin          = margin(10, 10, 10, 10)
  )
print(fig_deep)

ggsave(
  filename = "figure_prior_posterior_deep.pdf",
  plot = fig_deep,
  dpi = 300,
  units = "px",
  width = 3400,
  height = 2200
)
