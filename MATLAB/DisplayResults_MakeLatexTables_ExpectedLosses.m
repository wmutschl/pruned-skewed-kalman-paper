function DisplayResults_MakeLatexTables_ExpectedLosses(filename,ONLY_LATEX)
% DisplayResults_MakeLatexTables_ExpectedLosses(filename,ONLY_LATEX)
% -------------------------------------------------------------------------
% Collects and displays the results of the Monte-Carlo assessment of the 
% accuracy of filtered and smoothed states done in Accuracy_States.m
% Creates a log file with results and also the Latex code for the rows of
% Table 1 and Table 2 of the paper "Pruned Skewed Kalman Filter and Smoother:
% With Application to the Yield Curve" by Gaygysyz Guljanov, Willi Mutschler, Mark Trede
% -------------------------------------------------------------------------
% INPUTS
% - filename   [string]   load results created with Accuracy_States and saved in results/"filename".mat,
%                         for example: results_state_estimation_DGP1_T40_R2400_maca64
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

%% tables with expected losses
if ~ONLY_LATEX
    tblTitles.Loss = {'Gaussian'};
    for j=1:size(OPT.prune_tol,2)
        tblTitles.Loss = [tblTitles.Loss sprintf('CSN(%.d)',OPT.prune_tol(j))];
    end
end    
strKind = ["filt", "smooth"];
for j=1:length(strKind)
    if ~ONLY_LATEX
        fprintf('%s\nExpected Losses: %s\n',repmat('*',1,60),strKind(j))
    end
    MeanL = []; StdL = []; LowL = []; HighL = []; 
    RowNamesLoss = "";
    for strLoss = ["L1", "L2", "La"]
        if OPT.loss_fct.type.(strLoss)
            L = [RES.Loss.gauss.(strKind(j)).(strLoss) RES.Loss.csn.(strKind(j)).(strLoss)];
            MeanL = [MeanL; mean(L,1)];
            StdL = [StdL; std(L,[],1)];
            LowL  = [LowL; quantile(L,0.05,1)];
            HighL = [HighL; quantile(L,0.95,1)];
            RowNamesLoss = [RowNamesLoss sprintf("Mean %s",strLoss) sprintf("Std %s",strLoss) sprintf("Low %s",strLoss) sprintf("High %s",strLoss)];
        end
    end
    RowNamesLoss(1) = [];
    if ~ONLY_LATEX
        format long;
        disp(array2table([MeanL;StdL;LowL;HighL],'VariableNames',tblTitles.Loss,'RowNames',RowNamesLoss'));
        format short;
    end

    % Latex table
    fprintf('\n\nLatex Table Entries for %s:\n\n',strKind(j))
    if OPT.prune_tol(1) ~= 0
        MeanL = [MeanL(:,1) zeros(size(MeanL,1),1) MeanL(:,2:end)];
        LowL  = [LowL(:,1)  zeros(size(LowL,1),1)  LowL(:,2:end)];
        HighL = [HighL(:,1) zeros(size(HighL,1),1) HighL(:,2:end)];    
    end
    if strcmp(OPT.parameter_set,'DGP1')
        STRLOSS = ["L_1","L_2","L_a"];
        tbl_latex = [];
        for jj = 1:3
            strLoss = STRLOSS(jj);
            tbl_latex = [tbl_latex;
                sprintf('(1) & %d & $%s$ & $\\underset{[%.4f;%.4f]}{%.9f}$ & $\\underset{[%.4f;%.4f]}{%.9f}$ & $\\underset{[%.4f;%.4f]}{%.9f}$ & $\\underset{[%.4f;%.4f]}{%.9f}$ & $\\underset{[%.4f;%.4f]}{%.9f}$ \\\\'...
                ,OPT.sample_size...
                ,strLoss...
                ,LowL(jj,1), HighL(jj,1), MeanL(jj,1)...
                ,LowL(jj,2), HighL(jj,2), MeanL(jj,2)...
                ,LowL(jj,3), HighL(jj,3), MeanL(jj,3)...
                ,LowL(jj,4), HighL(jj,4), MeanL(jj,4)...
                ,LowL(jj,5), HighL(jj,5), MeanL(jj,5)...
                );];
        end
    elseif strcmp(OPT.parameter_set,'DGP2')
        strLoss = "L_2";
        jj=1;
        tbl_latex = [...
                sprintf('(2) & %d & $%s$ & $\\underset{[%.4f;%.4f]}{%.8f}$ & $\\underset{[%.4f;%.4f]}{%.8f}$ & $\\underset{[%.4f;%.4f]}{%.8f}$ & $\\underset{[%.4f;%.4f]}{%.8f}$ & $\\underset{[%.4f;%.4f]}{%.8f}$ \\\\'...
                ,OPT.sample_size...
                ,strLoss...
                ,LowL(jj,1), HighL(jj,1), MeanL(jj,1)...
                ,LowL(jj,2), HighL(jj,2), MeanL(jj,2)...
                ,LowL(jj,3), HighL(jj,3), MeanL(jj,3)...
                ,LowL(jj,4), HighL(jj,4), MeanL(jj,4)...
                ,LowL(jj,5), HighL(jj,5), MeanL(jj,5)...
                );];
    end
    disp(tbl_latex);
end
if ~ONLY_LATEX
    disp(['Total computing time : ' dynsec2hms(RES.time_Accuracy_States) ]);
end

diary off