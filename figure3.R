# Replicates Figure 3: "Probability density functions and cumulative distribution
# functions of a CSN distributed random variable with two skewness dimensions
# and the approximating pruned CSN distribution with one skewness dimenstions"
# =========================================================================
# Copyright © 2022-2023 Gaygysyz Guljanov, Willi Mutschler, Mark Trede
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
library(patchwork)
# set working directory
setwd(this.path::here())

source("R/csn_prune_params.R")

Sigma <- matrix(1,1,1)
Gamma <- matrix(c(6,0.1),2,1)
nu <- matrix(c(0,0),2,1)
Delta <- matrix(c(1,-0.1,-0.1,1),2,2)

P <- rbind(cbind(Sigma,Sigma%*%t(Gamma)),cbind(Gamma%*%Sigma, Delta+Gamma%*%Sigma%*%t(Gamma)))
R <- cov2cor(P)
R
PrunedParams <- pruned_csn_params(Sigma,Gamma,nu,Delta)

Sigma_pruned <- PrunedParams$Sigma
Gamma_pruned <- PrunedParams$Gamma
nu_pruned <- PrunedParams$nu
Delta_pruned <- PrunedParams$Delta

x <- seq(-1,4,length=200)
H <- data.frame(x=x,
                f1=dcsn(x,0,Sigma,        Gamma,        as.vector(nu),        Delta       ),
                f2=dcsn(x,0,Sigma_pruned, Gamma_pruned, as.vector(nu_pruned), Delta_pruned),
                c1=pcsn(x,0,Sigma,        Gamma,        as.vector(nu),        Delta       ),
                c2=pcsn(x,0,Sigma_pruned, Gamma_pruned, as.vector(nu_pruned), Delta_pruned)
                )
text_size = 12
title_size = 13

p1 <- ggplot(H, aes(x,f1)) +
        geom_line(linewidth=1.2)+
        geom_line(aes(y=f2),linetype="dashed",linewidth=1.2)+
        ggtitle("Probability distribution function")+
        ylab("")+
        theme(
          text = element_text(size=text_size),
          axis.title.y = element_text(size=text_size),
          axis.title.x = element_text(size=text_size),
          plot.title = element_text(size=title_size, hjust=0.5)
        )

p2 <- ggplot(H, aes(x,c1)) +
        geom_line(linewidth=1.2)+
        geom_line(aes(y=c2),linetype="dashed",linewidth=1.2)+
        ggtitle("Cumulative distribution function")+
        ylab("")+
        theme(
          text = element_text(size=text_size),
          axis.title.y = element_text(size=text_size),
          axis.title.x = element_text(size=text_size),
          plot.title = element_text(size=title_size, hjust=0.5)
        )

grid.arrange(p1,p2,ncol=2)
fig_3 <- (p1 + p2) + plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.box = "horizontal")
print(fig_3)
ggsave(filename = "plots/R/fig_3.pdf", plot = fig_3, dpi = 300, units = "px", width = 3000, height = 1000)
