function x = param_transform_bounded(y,lb,ub)
% transforms parameters y with unbounded support (used in optimization)
% into parameters x with bounded support (used in model)

% create indices
%idx_param_inf_inf = logical( isinf(lb) .*  isinf(ub) );
idx_param_lb_inf  = logical(~isinf(lb) .*  isinf(ub) );
idx_param_inf_ub  = logical( isinf(lb) .* ~isinf(ub) );
idx_param_lb_ub   = logical(~isinf(lb) .* ~isinf(ub) );

% for parameters x with support (-Inf Inf) do nothing, so we initialize x
x=y;

% for parameters x with support [lb Inf): x=exp(y)+lb (inverse transformation of y=log(x-lb) )
x(idx_param_lb_inf) = exp(y(idx_param_lb_inf)) + lb(idx_param_lb_inf) ;

% for parameters x with support (-Inf ub]: x=ub-exp(y) (inverse transformation of y=log(ub-x) )
x(idx_param_inf_ub) = ub(idx_param_inf_ub) - exp(y(idx_param_inf_ub));

% for parameters x with support (lb ub): (inverse transformation of y = logit((x-lb)/(ub-lb)) );
x(idx_param_lb_ub) = lb(idx_param_lb_ub) + (ub(idx_param_lb_ub)-lb(idx_param_lb_ub)).*inv_logit(y(idx_param_lb_ub));


function inv_logit_v = inv_logit(v)
% logistic sigmoid function defined for v in (-Inf,Inf); inverse of log-odds function
inv_logit_v = 1./(1+exp(v));
end

end % main function end