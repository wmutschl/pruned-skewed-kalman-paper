function Y = transform_to_unbounded(X,a,b)
% https://mc-stan.org/docs/reference-manual/logit-transform-jacobian.html
Y = logit((X-a)./(b-a));

function logit_u = logit(u)
    logit_u = log(u./(1-u));
end % logit

end % transform_to_unbounded