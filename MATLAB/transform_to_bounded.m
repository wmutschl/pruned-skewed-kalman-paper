function X = transform_to_bounded(Y,a,b)
% https://mc-stan.org/docs/reference-manual/logit-transform-jacobian.html
X = a + (b-a).*inv_logit(Y);

function u = inv_logit(v)
    u = 1./(1+exp(-v));
end %inv_logit

end % transform_to_unbounded