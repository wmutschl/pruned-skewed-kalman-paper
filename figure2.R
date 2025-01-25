# Replicates Figure 2: "Density functions of bivariate CSN distributions
# with different skewness parameters"
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
library(csn)
library(ggplot2)
library(gridExtra)
library(viridis)
library(latex2exp)
library(patchwork)
# set working directory
setwd(this.path::here())

source("R/csn_prune_params.R")

pltobj <- function(xx, mu, Sigma, Gamma, nu, Delta, ttl=""){
  text_size = 12
  title_size = 13
  z <- dcsn(xx, mu, Sigma, Gamma, nu, Delta)
  p <- ggplot(data.frame(x1=xx[,1],
                         x2=xx[,2],
                         f=z),
              aes(x1, x2, fill=f)) +
    geom_raster(interpolate=TRUE) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = TeX("$x_1$"), y = TeX("$x_2$")) +
    guides(fill="none") + 
    #scale_fill_gradientn(colours=viridis(10))+
    scale_fill_gradient(low = "white", high = "black") +
    ggtitle(ttl)+
    geom_vline(xintercept=0, colour="black")+
    geom_hline(yintercept=0, colour="black")+
    theme(
      text = element_text(size=text_size),
      axis.title.y = element_text(size=text_size),
      axis.title.x = element_text(size=text_size),
      plot.title = element_text(size=title_size, hjust=0.5)
    )
  return(p)
}

x1 <- seq(-1.2, 2, length=51)
x2 <- seq(-1.2, 2, length=53)
xx <- cbind(rep(x1,times=length(x2)), rep(x2, each=length(x1)))

mu_a = c(0,0); Sigma_a = matrix(c(1,0.7,0.7,1),2,2); Gamma_a = matrix(c(6,0,0, 6),2,2); nu_a = c( 0, 0); Delta_a = diag(2);
mu_b = c(0,0); Sigma_b = matrix(c(1,0.7,0.7,1),2,2); Gamma_b = matrix(c(6,0,0,-6),2,2); nu_b = c( 0, 0); Delta_b = diag(2);
mu_c = c(0,0); Sigma_c = matrix(c(1,0.7,0.7,1),2,2); Gamma_c = matrix(c(6,6,6, 6),2,2); nu_c = c( 0, 0); Delta_c = diag(2);
mu_d = c(0,0); Sigma_d = matrix(c(1,0.7,0.7,1),2,2); Gamma_d = matrix(c(6,0,0, 6),2,2); nu_d = c(-6,-6); Delta_d = diag(2);

p1 <- pltobj(xx, mu=mu_a, Sigma=Sigma_a, Gamma=Gamma_a, nu=nu_a, Delta=Delta_a, "(a)")
p2 <- pltobj(xx, mu=mu_b, Sigma=Sigma_b, Gamma=Gamma_b, nu=nu_b, Delta=Delta_b, "(b)")
p3 <- pltobj(xx, mu=mu_c, Sigma=Sigma_c, Gamma=Gamma_c, nu=nu_c, Delta=Delta_c, "(c)")
p4 <- pltobj(xx, mu=mu_d, Sigma=Sigma_d, Gamma=Gamma_d, nu=nu_d, Delta=Delta_d, "(d)")
grid.arrange(p1,p2,p3,p4,ncol=2)

fig_2_top    <- (p1 + p2) + plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.box = "horizontal")
fig_2_bottom <- (p3 + p4) + plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.box = "horizontal")
fig_2 <- (fig_2_top / fig_2_bottom)
print(fig_2)
ggsave(filename = "plots/R/fig_2.pdf", plot = fig_2, dpi = 300, units = "px", width = 3000, height = 2000)

