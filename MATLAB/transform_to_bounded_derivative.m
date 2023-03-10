function dX_dY = transform_to_bounded_derivative(Y,a,b)
% https://mc-stan.org/docs/reference-manual/logit-transform-jacobian.html
dX_dY = (b-a).*inv_logit(Y).*(1-inv_logit(Y));

function u = inv_logit(v)
    u = 1./(1+exp(-v));
end %inv_logit

end % transform_to_unbounded