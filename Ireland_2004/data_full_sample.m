load gpr.dat;
timeline=(1948.25:0.25:2003)';
%% Full sample
gobs = demean(gpr(:,1));
piobs = demean(gpr(:,2));
robs = demean(gpr(:,3));

ghat = demean(gpr(:,1));
pihat = demean(gpr(:,2));
rhat = demean(gpr(:,3));