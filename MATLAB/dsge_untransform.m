function [xparams_untransformed, bounds_untransformed] = dsge_untransform(xparams, bounds, estim_params_)
% function [xparams_untransformed, bounds_untransformed] = dsge_untransform(xparams, bounds, estim_params_)
% -------------------------------------------------------------------------
% create vectors with untransformed parameter values and untransformed bounds
% -------------------------------------------------------------------------
% INPUTS
% - xparams          [nparam x 1]   parameter vector with possibly transformed parameters
% - bounds           [nparam x 2]   lower and upper bounds (are infinite for transformed parameters)
% - estim_params_    [structure]    information on estimated parameters
% -------------------------------------------------------------------------
% OUTPUTS
% - xparams_untransformed   [nparam x 1]   parameter vector with untransformed (original) parameters
% - bounds_untransformed    [nparam x 2]   lower and upper bounds for untransformed (original) parameters
% =========================================================================
% Copyright (C) 2023 Willi Mutschler
%
% This is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% -------------------------------------------------------------------------
% This file is part of the replication files for the paper "Pruned Skewed
% Kalman Filter and Smoother: Pruned Skewed Kalman Filter and Smoother:
% With Applications to the Yield Curve and Asymmetric Monetary Policy Shocks"
% by Gaygysyz Guljanov, Willi Mutschler, Mark Trede
% =========================================================================
xparams_untransformed = xparams;
bounds_untransformed = bounds;
for j=1:size(xparams,2)
    for jp=1:estim_params_.nsx
        if estim_params_.transformed(jp)
            idx = jp;
            if j==1
                bounds_untransformed(idx,1) = estim_params_.skew_exo(jp,3); bounds_untransformed(idx,2) = estim_params_.skew_exo(jp,4);
            end
            xparams_untransformed(idx,j) = transform_to_bounded(xparams(idx,j),bounds_untransformed(idx,1),bounds_untransformed(idx,2));
        end
    end
    for jp=1:estim_params_.nvx
        if estim_params_.transformed(estim_params_.nsx+jp)
            idx = estim_params_.nsx + jp;
            if j==1
                bounds_untransformed(idx,1) = estim_params_.var_exo(jp,3); bounds_untransformed(idx,2) = estim_params_.var_exo(jp,4);
            end
            xparams_untransformed(idx,j) = transform_to_bounded(xparams(idx,j),bounds_untransformed(idx,1),bounds_untransformed(idx,2));
        end
    end
    for jp=1:estim_params_.nvn
        if estim_params_.transformed(estim_params_.nsx+estim_params_.nvx+jp)
            idx = estim_params_.nsx + estim_params_.nvx + jp;
            if j==1
                bounds_untransformed(idx,1) = estim_params_.var_endo(jp,3); bounds_untransformed(idx,2) = estim_params_.var_endo(jp,4);
            end
            xparams_untransformed(idx,j) = transform_to_bounded(xparams(idx,j),bounds_untransformed(idx,1),bounds_untransformed(idx,2));
        end
    end
    for jp=1:estim_params_.np
        if estim_params_.transformed(estim_params_.nsx+estim_params_.nvx+estim_params_.nvn+jp)
            idx = estim_params_.nsx + estim_params_.nvx + estim_params_.nvn + jp;
            if j==1
                bounds_untransformed(idx,1) = estim_params_.param_vals(jp,3); bounds_untransformed(idx,2) = estim_params_.param_vals(jp,4);
            end
            xparams_untransformed(idx,j) = transform_to_bounded(xparams(idx,j),bounds_untransformed(idx,1),bounds_untransformed(idx,2));
        end
    end
end