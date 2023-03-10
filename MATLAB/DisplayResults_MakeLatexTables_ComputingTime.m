function DisplayResults_MakeLatexTables_ComputingTime(filename,ONLY_LATEX)
% DisplayResults_MakeLatexTables_ComputingTime(filename,ONLY_LATEX)
% -------------------------------------------------------------------------
% Collects and displays the results of the Monte-Carlo assessment of the 
% computing time of the log-likelihood function with both the Gaussian
% as well as Pruned Skewed Kalman filter in Computing_Time_Log_Likliheood.m
% Creates a log file with results and also the Latex code for the rows of
% Table 3 of the paper "Pruned Skewed Kalman Filter and Smoother: With Application
% to the Yield Curve" by Gaygysyz Guljanov, Willi Mutschler, Mark Trede
% -------------------------------------------------------------------------
% INPUTS
% - filename   [string]   load results created with Computing_Time_Log_Likliheood.m and saved in results/"filename".mat,
%                         for example: results_computing_time_DGP1_T50_R1000_maca64
% - ONLY_LATEX [boolean]  only print the latex rows in the command window
% -------------------------------------------------------------------------
% OUTPUTS
% - results are displayed in the command window
% - latex codes are displayed in the command window
% - a log file "filename.log" is created in the results folder
% =========================================================================
% Copyright (C) 2022-2023 Willi Mutschler
%
% This is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% -------------------------------------------------------------------------
% This file is part of the replication files for the paper "Pruned Skewed
% Kalman Filter and Smoother: With Application to the Yield Curve" by
% Gaygysyz Guljanov, Willi Mutschler, Mark Trede (2022)
% =========================================================================
if nargin < 2
    ONLY_LATEX = false;
end

load(['results/' filename],'OPT','PARAMS','RES');

%% start log file
if ~ONLY_LATEX
    diary off
    OPT.logfile = ['results/' filename,'.log'];
    if exist(OPT.logfile, 'file')
        delete(OPT.logfile)
    end
    diary(OPT.logfile)
end

if ~ONLY_LATEX
    format long
    disp(repmat('*',1,60));
    % print settings to logfile
    disp(OPT)
    disp(repmat('*',1,60));
    disp('Parameter Matrix G:');              disp(PARAMS.G);
    disp('Parameter Matrix F:');              disp(PARAMS.F);
    disp('Parameter Matrix R:');              disp(PARAMS.R);
    disp('Parameter Matrix mu_eps:');         disp(PARAMS.mu_eps);
    disp('Parameter Matrix Sigma_eps:');      disp(PARAMS.Sigma_eps);
    disp('Parameter Matrix mu_eta:');         disp(PARAMS.mu_eta);
    disp('Parameter Matrix Sigma_eta:');      disp(PARAMS.Sigma_eta);
    disp('Parameter Matrix nu_eta:');         disp(PARAMS.nu_eta);
    disp('Parameter Matrix Gamma_eta:');      disp(PARAMS.Gamma_eta);
    disp('Parameter Matrix Delta_eta:');      disp(PARAMS.Delta_eta);
    disp('Parameter Matrix lambda_eta:');     disp(PARAMS.lambda_eta);
    disp('Parameter Matrix sqrt_Sigma_eta:'); disp(PARAMS.sqrt_Sigma_eta);
    disp(repmat('*',1,60));
    format short
end

%% table with summary statistics
if ~ONLY_LATEX
    tblTitles.SampleStats = {};
    for j=1:size(PARAMS.F,1)
        tblTitles.SampleStats = [tblTitles.SampleStats sprintf('mean(y%d)',j) sprintf('sd(y%d)',j) sprintf('skew(y%d)',j)];
    end
    fprintf('%s\n%sSUMMARY STATISTICS\n',repmat('*',1,60),repmat(' ',1,10));
    disp(array2table(RES.SampleStats,'VariableNames',tblTitles.SampleStats));
end

%% table with computing time
Time     = [RES.Time.gauss RES.Time.csn]*1000; % in miliseconds
MeanTime = mean(Time,1);
LowTime  = quantile(Time,0.05,1);
HighTime = quantile(Time,0.95,1);
RowNamesLoss = [sprintf("Mean") sprintf("Low") sprintf("High")];

if ~ONLY_LATEX
    fprintf('%s\nTime to compute the log-likelihood once (in ms)\n',repmat('*',1,60))
    tblTitles.Time = {'Gaussian'};
    for j=1:size(OPT.prune_tol,2)
        tblTitles.Time = [tblTitles.Time sprintf('CSN(%.d)',OPT.prune_tol(j))];
    end
    format long;
    disp(array2table([MeanTime;LowTime;HighTime],'VariableNames',tblTitles.Time,'RowNames',RowNamesLoss'));
    format short;
    disp(['Total computing time : ' dynsec2hms(RES.time_Computing_Time) ]);
end

%% Latex table
if ~ONLY_LATEX
    fprintf('\n\nLatex Table Entries:\n\n')
end
if OPT.prune_tol(1) ~= 0
    MeanTime = [MeanTime(:,1) zeros(size(MeanTime,1),1) MeanTime(:,2:end)];
    LowTime  = [LowTime(:,1)  zeros(size(LowTime,1),1)  LowTime(:,2:end)];
    HighTime = [HighTime(:,1) zeros(size(HighTime,1),1) HighTime(:,2:end)];
end
tbl_latex = [...
            sprintf('%s & %d & $\\underset{[%.2f;%.2f]}{%.4f}$ & $\\underset{[%.2f;%.2f]}{%.4f}$ & $\\underset{[%.2f;%.2f]}{%.4f}$ & $\\underset{[%.2f;%.2f]}{%.4f}$ & $\\underset{[%.2f;%.2f]}{%.4f}$ \\\\'...
            ,OPT.parameter_set...
            ,OPT.sample_size...
            ,LowTime(1,1), HighTime(1,1), MeanTime(1,1)...
            ,LowTime(1,2), HighTime(1,2), MeanTime(1,2)...
            ,LowTime(1,3), HighTime(1,3), MeanTime(1,3)...
            ,LowTime(1,4), HighTime(1,4), MeanTime(1,4)...
            ,LowTime(1,5), HighTime(1,5), MeanTime(1,5)...
            );];
disp(tbl_latex);

if ~ONLY_LATEX
    diary off
end