function y = param_transform_unbounded(x,lb,ub)
% transforms parameters x with bounded support (used in model) into
% parameters y with unbounded support (used in optimization)

% create indices
%idx_param_inf_inf = logical( isinf(lb) .*  isinf(ub) );
idx_param_lb_inf  = logical(~isinf(lb) .*  isinf(ub) );
idx_param_inf_ub  = logical( isinf(lb) .* ~isinf(ub) );
idx_param_lb_ub   = logical(~isinf(lb) .* ~isinf(ub) );

% for parameters x with support (-Inf Inf) do nothing, so we initialize y
y=x;

% for parameters x with support [lb Inf): y = log(x-lb)
y(idx_param_lb_inf) = log(x(idx_param_lb_inf) - lb(idx_param_lb_inf));

% for parameters x with support (-Inf ub]: y = log(ub-x)
y(idx_param_inf_ub) = log(ub(idx_param_inf_ub)-x(idx_param_inf_ub));

% for parameters x with support (lb ub): y = logit((x-lb)/(ub-lb));
y(idx_param_lb_ub) = logit( (x(idx_param_lb_ub)-lb(idx_param_lb_ub))./(ub(idx_param_lb_ub)-lb(idx_param_lb_ub)) );


function logit_u = logit(u)
% log-odds function for u in (0,1); inverse of logistic sigmoid function
logit_u = log(u./(1-u));
end

end%main function end