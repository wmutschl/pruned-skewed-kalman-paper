function xstderr = standard_errors_inverse_hessian(objfct, xparams_untransformed, bounds_untransformed, varargin)
xstderr = nan(size(xparams_untransformed));
parfor jx = 1:size(xparams_untransformed,2)
    xparams = xparams_untransformed(:,jx);
    gstep = max(abs(xparams), sqrt(1e-3)*ones(length(xparams), 1))*eps^(1/6);
    H = get_hessian_with_bounds(objfct, xparams, gstep, bounds_untransformed,    varargin{:});
    V = inv(reshape(H,length(xparams),length(xparams)));
    xstderr(:,jx) = sqrt(diag(V));
end