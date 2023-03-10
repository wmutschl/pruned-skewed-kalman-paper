function SE = get_standard_errors(objfct,xparams,bounds,datamat,estim_params_,options_,M_)
    H = reshape(get_hessian(objfct,xparams,[1e-3;1.0],bounds,datamat,estim_params_,options_,M_),estim_params_.ntot,estim_params_.ntot);
    V = inv(H);
    SE = sqrt(diag(V));

for jp=1:estim_params0.ntot
    if estim_params0.transformed(jp)
        dX_dY = transform_to_bounded_derivative(xparams0Final(11),0,1)
        (dX_dY*SE_stage0(jp))^2
end
    