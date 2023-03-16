function [xparams_untransformed, bounds_untransformed] = dsge_untransform(xparams, bounds, estim_params_)
% create vectors with untransform parameter values and bounds
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