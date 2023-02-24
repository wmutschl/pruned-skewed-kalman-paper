function [STEADY_STATE,PARAM,error_indicator] = ireland_2004_steady_state(MODEL,PARAM)
error_indicator = 0; % initialize no error

% read-out parameters

% compute steady-state
a = 0;
e = 0;
z = 0;
x = 0;
pihat = 0;
yhat = 0;
ghat = 0;
rhat = 0;

eta_a = 0;
eta_e = 0;
eta_z = 0;
eta_r = 0;

% write to output structure, important keep same ordering as in MODEL.endo_names and MODEL.exo_names
for j = 1:MODEL.endo_nbr
    varname = MODEL.endo_names(j,1);
    STEADY_STATE.(varname) = eval(varname);
end

for j = 1:MODEL.exo_nbr
    varname = MODEL.exo_names(j,1);
    STEADY_STATE.(varname) = eval(varname);
end


% check if anything went wrong
if any(isnan(struct2values(STEADY_STATE))) || any(isinf(struct2values(STEADY_STATE))) || any(imag(struct2values(STEADY_STATE)))
    error_indicator = 1;    
    warning('Something wrong with the steady-state computations due to nan, inf or imaginary numbers');
    return
end
% check residuals
if any( feval(str2func(MODEL.name + '_f'),STEADY_STATE,PARAM) > 1e-7 )
    error_indicator = 1;
    warning('Something wrong with the steady-state computations as the residuals are nonzero');
    return
end

end% main function end