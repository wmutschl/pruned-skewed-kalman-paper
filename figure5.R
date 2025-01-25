# Replicates Figure 5: ""
# =========================================================================
# Copyright © 2025 Gaygysyz Guljanov, Willi Mutschler, Mark Trede
#
# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# -------------------------------------------------------------------------
# This file is part of the replication files for the paper "Pruned Skewed
# Kalman Filter and Smoother: Pruned Skewed Kalman Filter and Smoother:
# With Applications to the Yield Curve and Asymmetric Monetary Policy Shocks"
# by Gaygysyz Guljanov, Willi Mutschler, Mark Trede
# =========================================================================
library(ggplot2)
library(gridExtra)
library(csn)
library(latex2exp)
library(R.matlab)
library(dplyr)
library(tidyr)
# set working directory
setwd(this.path::here())

source("R/csn_prune_params.R")

irfs <- readMat("ireland2004_irfs/Output/ireland2004_irfs.mat")
gaussian_neg <- irfs_gaussian_neg[, , 1];
gaussian_pos <- irfs_gaussian_pos[, , 1]
csn_neg      <- irfs_csn_neg[, , 1]
csn_pos      <- irfs_csn_pos[, , 1]

# monetary policy shocks
irfLength = 15
df <- bind_rows(

  # interest rate & preference shock
  tibble(time = 1:irfLength, distribution = "Gaussian", sign = "negative", shock = "preference", variable = "interest rate", value = as.numeric(gaussian_neg$rhat.eta.a)),
  tibble(time = 1:irfLength, distribution = "Gaussian", sign = "positive", shock = "preference", variable = "interest rate", value = as.numeric(gaussian_pos$rhat.eta.a)),
  tibble(time = 1:irfLength, distribution = "CSN",      sign = "negative", shock = "preference", variable = "interest rate", value = as.numeric(csn_neg$rhat.eta.a)     ),
  tibble(time = 1:irfLength, distribution = "CSN",      sign = "positive", shock = "preference", variable = "interest rate", value = as.numeric(csn_pos$rhat.eta.a)     ),
  # interest rate & cost-push shock
  tibble(time = 1:irfLength, distribution = "Gaussian", sign = "negative", shock = "cost-push", variable = "interest rate", value = as.numeric(gaussian_neg$rhat.eta.e)),
  tibble(time = 1:irfLength, distribution = "Gaussian", sign = "positive", shock = "cost-push", variable = "interest rate", value = as.numeric(gaussian_pos$rhat.eta.e)),
  tibble(time = 1:irfLength, distribution = "CSN",      sign = "negative", shock = "cost-push", variable = "interest rate", value = as.numeric(csn_neg$rhat.eta.e)     ),
  tibble(time = 1:irfLength, distribution = "CSN",      sign = "positive", shock = "cost-push", variable = "interest rate", value = as.numeric(csn_pos$rhat.eta.e)     ),
  # interest rate & productivity shock
  tibble(time = 1:irfLength, distribution = "Gaussian", sign = "negative", shock = "productivity", variable = "interest rate", value = as.numeric(gaussian_neg$rhat.eta.z)),
  tibble(time = 1:irfLength, distribution = "Gaussian", sign = "positive", shock = "productivity", variable = "interest rate", value = as.numeric(gaussian_pos$rhat.eta.z)),
  tibble(time = 1:irfLength, distribution = "CSN",      sign = "negative", shock = "productivity", variable = "interest rate", value = as.numeric(csn_neg$rhat.eta.z)     ),
  tibble(time = 1:irfLength, distribution = "CSN",      sign = "positive", shock = "productivity", variable = "interest rate", value = as.numeric(csn_pos$rhat.eta.z)     ),
  # interest rate & monetary policy shock
  tibble(time = 1:irfLength, distribution = "Gaussian", sign = "negative", shock = "monetary policy", variable = "interest rate", value = as.numeric(gaussian_neg$rhat.eta.r)),
  tibble(time = 1:irfLength, distribution = "Gaussian", sign = "positive", shock = "monetary policy", variable = "interest rate", value = as.numeric(gaussian_pos$rhat.eta.r)),
  tibble(time = 1:irfLength, distribution = "CSN",      sign = "negative", shock = "monetary policy", variable = "interest rate", value = as.numeric(csn_neg$rhat.eta.r)     ),
  tibble(time = 1:irfLength, distribution = "CSN",      sign = "positive", shock = "monetary policy", variable = "interest rate", value = as.numeric(csn_pos$rhat.eta.r)     ),
  
  # output gap & preference shock
  tibble(time = 1:irfLength, distribution = "Gaussian", sign = "negative", shock = "preference", variable = "output gap", value = as.numeric(gaussian_neg$xhat.eta.a)),
  tibble(time = 1:irfLength, distribution = "Gaussian", sign = "positive", shock = "preference", variable = "output gap", value = as.numeric(gaussian_pos$xhat.eta.a)),
  tibble(time = 1:irfLength, distribution = "CSN",      sign = "negative", shock = "preference", variable = "output gap", value = as.numeric(csn_neg$xhat.eta.a)     ),
  tibble(time = 1:irfLength, distribution = "CSN",      sign = "positive", shock = "preference", variable = "output gap", value = as.numeric(csn_pos$xhat.eta.a)     ),
  # output gap & cost-push shock
  tibble(time = 1:irfLength, distribution = "Gaussian", sign = "negative", shock = "cost-push", variable = "output gap", value = as.numeric(gaussian_neg$xhat.eta.e)),
  tibble(time = 1:irfLength, distribution = "Gaussian", sign = "positive", shock = "cost-push", variable = "output gap", value = as.numeric(gaussian_pos$xhat.eta.e)),
  tibble(time = 1:irfLength, distribution = "CSN",      sign = "negative", shock = "cost-push", variable = "output gap", value = as.numeric(csn_neg$xhat.eta.e)     ),
  tibble(time = 1:irfLength, distribution = "CSN",      sign = "positive", shock = "cost-push", variable = "output gap", value = as.numeric(csn_pos$xhat.eta.e)     ),
  # output gap & productivity shock
  tibble(time = 1:irfLength, distribution = "Gaussian", sign = "negative", shock = "productivity", variable = "output gap", value = as.numeric(gaussian_neg$xhat.eta.z)),
  tibble(time = 1:irfLength, distribution = "Gaussian", sign = "positive", shock = "productivity", variable = "output gap", value = as.numeric(gaussian_pos$xhat.eta.z)),
  tibble(time = 1:irfLength, distribution = "CSN",      sign = "negative", shock = "productivity", variable = "output gap", value = as.numeric(csn_neg$xhat.eta.z)     ),
  tibble(time = 1:irfLength, distribution = "CSN",      sign = "positive", shock = "productivity", variable = "output gap", value = as.numeric(csn_pos$xhat.eta.z)     ),
  # output gap & monetary policy shock
  tibble(time = 1:irfLength, distribution = "Gaussian", sign = "negative", shock = "monetary policy", variable = "output gap", value = as.numeric(gaussian_neg$xhat.eta.r)),
  tibble(time = 1:irfLength, distribution = "Gaussian", sign = "positive", shock = "monetary policy", variable = "output gap", value = as.numeric(gaussian_pos$xhat.eta.r)),
  tibble(time = 1:irfLength, distribution = "CSN",      sign = "negative", shock = "monetary policy", variable = "output gap", value = as.numeric(csn_neg$xhat.eta.r)     ),
  tibble(time = 1:irfLength, distribution = "CSN",      sign = "positive", shock = "monetary policy", variable = "output gap", value = as.numeric(csn_pos$xhat.eta.r)     ),
  
  # output growth & preference shock
  tibble(time = 1:irfLength, distribution = "Gaussian", sign = "negative", shock = "preference", variable = "output growth", value = as.numeric(gaussian_neg$ghat.eta.a)),
  tibble(time = 1:irfLength, distribution = "Gaussian", sign = "positive", shock = "preference", variable = "output growth", value = as.numeric(gaussian_pos$ghat.eta.a)),
  tibble(time = 1:irfLength, distribution = "CSN",      sign = "negative", shock = "preference", variable = "output growth", value = as.numeric(csn_neg$ghat.eta.a)     ),
  tibble(time = 1:irfLength, distribution = "CSN",      sign = "positive", shock = "preference", variable = "output growth", value = as.numeric(csn_pos$ghat.eta.a)     ),
  # output growth & cost-push shock
  tibble(time = 1:irfLength, distribution = "Gaussian", sign = "negative", shock = "cost-push", variable = "output growth", value = as.numeric(gaussian_neg$ghat.eta.e)),
  tibble(time = 1:irfLength, distribution = "Gaussian", sign = "positive", shock = "cost-push", variable = "output growth", value = as.numeric(gaussian_pos$ghat.eta.e)),
  tibble(time = 1:irfLength, distribution = "CSN",      sign = "negative", shock = "cost-push", variable = "output growth", value = as.numeric(csn_neg$ghat.eta.e)     ),
  tibble(time = 1:irfLength, distribution = "CSN",      sign = "positive", shock = "cost-push", variable = "output growth", value = as.numeric(csn_pos$ghat.eta.e)     ),
  # output growth & productivity shock
  tibble(time = 1:irfLength, distribution = "Gaussian", sign = "negative", shock = "productivity", variable = "output growth", value = as.numeric(gaussian_neg$ghat.eta.z)),
  tibble(time = 1:irfLength, distribution = "Gaussian", sign = "positive", shock = "productivity", variable = "output growth", value = as.numeric(gaussian_pos$ghat.eta.z)),
  tibble(time = 1:irfLength, distribution = "CSN",      sign = "negative", shock = "productivity", variable = "output growth", value = as.numeric(csn_neg$ghat.eta.z)     ),
  tibble(time = 1:irfLength, distribution = "CSN",      sign = "positive", shock = "productivity", variable = "output growth", value = as.numeric(csn_pos$ghat.eta.z)     ),
  # output growth & monetary policy shock
  tibble(time = 1:irfLength, distribution = "Gaussian", sign = "negative", shock = "monetary policy", variable = "output growth", value = as.numeric(gaussian_neg$ghat.eta.r)),
  tibble(time = 1:irfLength, distribution = "Gaussian", sign = "positive", shock = "monetary policy", variable = "output growth", value = as.numeric(gaussian_pos$ghat.eta.r)),
  tibble(time = 1:irfLength, distribution = "CSN",      sign = "negative", shock = "monetary policy", variable = "output growth", value = as.numeric(csn_neg$ghat.eta.r)     ),
  tibble(time = 1:irfLength, distribution = "CSN",      sign = "positive", shock = "monetary policy", variable = "output growth", value = as.numeric(csn_pos$ghat.eta.r)     ),
  
  # inflation & preference shock
  tibble(time = 1:irfLength, distribution = "Gaussian", sign = "negative", shock = "preference", variable = "inflation", value = as.numeric(gaussian_neg$pihat.eta.a)),
  tibble(time = 1:irfLength, distribution = "Gaussian", sign = "positive", shock = "preference", variable = "inflation", value = as.numeric(gaussian_pos$pihat.eta.a)),
  tibble(time = 1:irfLength, distribution = "CSN",      sign = "negative", shock = "preference", variable = "inflation", value = as.numeric(csn_neg$pihat.eta.a)     ),
  tibble(time = 1:irfLength, distribution = "CSN",      sign = "positive", shock = "preference", variable = "inflation", value = as.numeric(csn_pos$pihat.eta.a)     ),
  # inflation & productivity shock
  tibble(time = 1:irfLength, distribution = "Gaussian", sign = "negative", shock = "cost-push", variable = "inflation", value = as.numeric(gaussian_neg$pihat.eta.e)),
  tibble(time = 1:irfLength, distribution = "Gaussian", sign = "positive", shock = "cost-push", variable = "inflation", value = as.numeric(gaussian_pos$pihat.eta.e)),
  tibble(time = 1:irfLength, distribution = "CSN",      sign = "negative", shock = "cost-push", variable = "inflation", value = as.numeric(csn_neg$pihat.eta.e)     ),
  tibble(time = 1:irfLength, distribution = "CSN",      sign = "positive", shock = "cost-push", variable = "inflation", value = as.numeric(csn_pos$pihat.eta.e)     ),
  # inflation & productivity shock
  tibble(time = 1:irfLength, distribution = "Gaussian", sign = "negative", shock = "productivity", variable = "inflation", value = as.numeric(gaussian_neg$pihat.eta.z)),
  tibble(time = 1:irfLength, distribution = "Gaussian", sign = "positive", shock = "productivity", variable = "inflation", value = as.numeric(gaussian_pos$pihat.eta.z)),
  tibble(time = 1:irfLength, distribution = "CSN",      sign = "negative", shock = "productivity", variable = "inflation", value = as.numeric(csn_neg$pihat.eta.z)     ),
  tibble(time = 1:irfLength, distribution = "CSN",      sign = "positive", shock = "productivity", variable = "inflation", value = as.numeric(csn_pos$pihat.eta.z)     ),
  # inflation & monetary policy shock
  tibble(time = 1:irfLength, distribution = "Gaussian", sign = "negative", shock = "monetary policy", variable = "inflation", value = as.numeric(gaussian_neg$pihat.eta.r)),
  tibble(time = 1:irfLength, distribution = "Gaussian", sign = "positive", shock = "monetary policy", variable = "inflation", value = as.numeric(gaussian_pos$pihat.eta.r)),
  tibble(time = 1:irfLength, distribution = "CSN",      sign = "negative", shock = "monetary policy", variable = "inflation", value = as.numeric(csn_neg$pihat.eta.r)     ),
  tibble(time = 1:irfLength, distribution = "CSN",      sign = "positive", shock = "monetary policy", variable = "inflation", value = as.numeric(csn_pos$pihat.eta.r)     ),
  
)

ggplot(df, aes(x = time, y = value, color = sign, linetype = distribution)) +
  geom_line(size = 1) +
  facet_grid(variable ~ shock, scales = "free_y") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Impulse Response Functions (IRFs)",
    x = "Time Horizon",
    y = "Response",
    color = "Shock Sign",
    linetype = "Distribution"
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "bottom"
  ) +
  scale_color_manual(values = c("negative" = "red", "positive" = "blue")) +
  scale_linetype_manual(values = c("Gaussian" = "solid", "CSN" = "dashed"))

df_monetary <- df %>%
  filter(shock == "monetary policy")

# Create the IRF plot
ggplot(df_monetary, aes(x = time, y = value, color = sign, linetype = distribution)) +
  geom_line(size = 1) +
  facet_wrap(~ variable, ncol = 1, scales = "free_y") +  # One panel per variable
  theme_minimal(base_size = 14) +
  labs(
    title = "Impulse Response Functions (IRFs) to Monetary Policy Shocks",
    subtitle = "Comparing Negative and Positive Shocks",
    x = "Time Horizon",
    y = "Response",
    color = "Shock Sign",
    linetype = "Distribution"
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12)
  ) +
  scale_color_manual(values = c("negative" = "red", "positive" = "blue")) +
  scale_linetype_manual(values = c("Gaussian" = "solid", "CSN" = "dashed"))
