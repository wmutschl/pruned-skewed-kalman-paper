load data/ireland2004_gpr.dat;
DATA_SAMPLE = "POST_1980";
timeline = (1948.25:0.25:2003)';
if DATA_SAMPLE == "PRE_1980"
    ghat  = demean(ireland2004_gpr(timeline<1980,1));
    pihat = demean(ireland2004_gpr(timeline<1980,2));
    rhat  = demean(ireland2004_gpr(timeline<1980,3));
elseif DATA_SAMPLE == "POST_1980"
    ghat  = demean(ireland2004_gpr(timeline>=1980,1));
    pihat = demean(ireland2004_gpr(timeline>=1980,2));
    rhat  = demean(ireland2004_gpr(timeline>=1980,3));
else % FULL SAMPLE
    ghat  = demean(ireland2004_gpr(:,1));
    pihat = demean(ireland2004_gpr(:,2));
    rhat  = demean(ireland2004_gpr(:,3));
end