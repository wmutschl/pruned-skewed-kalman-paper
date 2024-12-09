# Replicates Figure 1: "Density functions of univariate CSN distributions with
# different skewness parameters"
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
source("R/csn_prune_params.R")

mu <- 0
nu <- 0
Sigma <- matrix(1,1,1)
Delta <- matrix(1,1,1)

x <- seq(-4,4,length=200)
f1 <- dcsn(matrix(x,200,1),mu,Sigma,gamma=matrix(-5,1,1),nu,Delta)
f2 <- dcsn(matrix(x,200,1),mu,Sigma,gamma=matrix(0,1,1),nu,Delta)
f3 <- dcsn(matrix(x,200,1),mu,Sigma,gamma=matrix(5,1,1),nu,Delta)
f4 <- dcsn(matrix(x,200,1),mu,Sigma,gamma=matrix(5,1,1),-8,Delta)

p1 <- ggplot(data.frame(x=x,f=f1),aes(x,f))+
    geom_line()+
    ggtitle(TeX("(a) \\Gamma=-5, \\nu=0"))+
    theme(text = element_text(size=30), plot.title = element_text(size = 40))+
    ylab("density")
p2 <- ggplot(data.frame(x=x,f=f2),aes(x,f))+
    geom_line()+
    ggtitle(TeX("(b) \\Gamma=0, \\nu=0"))+
    theme(text = element_text(size=30), plot.title = element_text(size = 40))+
    ylab("")
p3 <- ggplot(data.frame(x=x,f=f3),aes(x,f))+
    geom_line()+
    ggtitle(TeX("(c) \\Gamma=5, \\nu=0"))+
    theme(text = element_text(size=30), plot.title = element_text(size = 40))+
    ylab("")
p4 <- ggplot(data.frame(x=x,f=f4),aes(x,f))+
    geom_line()+
    ggtitle(TeX("(d) \\Gamma=5, \\nu=-8"))+
    theme(text = element_text(size=30), plot.title = element_text(size = 40))+  
    ylab("")

grid.arrange(p1, p2, p3, p4, ncol=2)


# Bivariate CSN density examples

pltobj <- function(xx, mu, Sigma, D, nu, Delta, ttl=""){
    z <- dcsn(xx, mu, Sigma, D, nu, Delta)
    p <- ggplot(data.frame(x1=xx[,1],
                           x2=xx[,2],
                           f=z),
                aes(x1, x2, fill=f)) +
        geom_raster(interpolate=TRUE)+
        scale_x_continuous("X1", expand = c(0, 0)) +
        scale_y_continuous("X2", expand = c(0, 0)) +
        guides(fill="none") + 
        scale_fill_gradientn(colours=viridis(10))+
        ggtitle(ttl)+
        geom_vline(xintercept=0, colour="grey")+
        geom_hline(yintercept=0, colour="grey")
    return(p)
}

x1 <- seq(-1.2, 2, length=51)
x2 <- seq(-1.2, 2, length=53)
xx <- cbind(rep(x1,times=length(x2)), rep(x2, each=length(x1)))

p1 <- pltobj(xx,
             mu=c(0,0), 
             Sigma=matrix(c(1,0.7,0.7,1),2,2), 
             D=matrix(c(6,0,0,6),2,2), 
             nu=c(0,0), 
             Delta=diag(2),
             "(a)")

p2 <- pltobj(xx,
             mu=c(0,0), 
             Sigma=matrix(c(1,0.7,0.7,1),2,2), 
             D=matrix(c(6,0,0,-6),2,2), 
             nu=c(0,0), 
             Delta=diag(2),
             "(b)")

p3 <- pltobj(xx,
             mu=c(0,0), 
             Sigma=matrix(c(1,0.7,0.7,1),2,2), 
             D=matrix(c(6,6,6,6),2,2), 
             nu=c(0,0), 
             Delta=diag(2),
             "(c)")

p4 <- pltobj(xx,
             mu=c(0,0), 
             Sigma=matrix(c(1,0.7,0.7,1),2,2), 
             D=matrix(c(6,0,0,6),2,2), 
             nu=c(-6,-6), 
             Delta=diag(2),
             "(d)")

opendev("linux",6,6)
grid.arrange(p1,p2,p3,p4,ncol=2)
savedev("linux",file="~/Schreibtisch/biv_pdf.pdf")
