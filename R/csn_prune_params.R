pruned_csn_params <- function(Sigma,Gamma,nu,Delta, tol=0.1){
    # pruned_csn_params <- function(Sigma,Gamma,nu,Delta, tol=0.1)
    # -------------------------------------------------------------------------
    # compute pruned parameters of CSN(mu,Sigma,Gamma,nu,Delta) distribution
    # =========================================================================
    # INPUTS
    #   - Sigma: scale parameter of CSN distribution, [p x p] matrix
    #   - Gamma: first skewness parameter of CSN distribution, [q x p] matrix
    #   - nu   : second skewness parameter of CSN distribution, [q x 1] vector
    #   - Delta: third skewness parameter of CSN distribution, [q x 1] matrix
    #   - tol  : threshold at which correlations are pruned, scalar between 0 and 1
    # -------------------------------------------------------------------------
    # OUTPUTS
    #   - pruned_csn_params: list of pruned parameters Sigma, Gamma, nu, Delta
    # =========================================================================
    # Copyright © 2022-2023 Gaygysyz Guljanov, Willi Mutschler
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
    p <- dim(Sigma)[1]
    q <- dim(Delta)[1]
    P1 <- Sigma
    P2 <- Sigma %*% t(Gamma)
    P3 <- Gamma %*% Sigma
    P4 <- Delta+Gamma %*% Sigma %*% t(Gamma)
    P <- rbind(cbind(P1,P2), cbind(P3,P4))
    maxval <- apply(cov2cor(P)[(p+1):(p+q), 1:p, drop=F], 1, max)
    cutvector <- (maxval < tol)
    if(all(cutvector)){
        Gamma <- matrix(0,1,1)
        nu <- 0
        Delta <- matrix(1,1,1)
    } else if(any(cutvector)){
        P3 <- P3[!cutvector,,drop=F]
        P4 <- P4[!cutvector,!cutvector,drop=F]
        Gamma <- P3 %*% solve(Sigma)
        nu <- nu[!cutvector]
        Delta <- P4 - P3 %*% t(Gamma)
    }
    return(list(Sigma=Sigma,
                Gamma=Gamma,
                nu=nu,
                Delta=Delta))
}