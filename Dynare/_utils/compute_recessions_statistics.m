function [meanDuration,meanOutputLoss,recessionNumber,T] = compute_recessions_statistics(ghat, cutoff, minQuarters)
% function [meanDuration,meanOutputLoss,recessionNumber,T] = compute_recessions_statistics(ghat, cutoff, minQuarters)
% -------------------------------------------------------------------------
% Computes statistics on recessions based on simulated output growth rate ghat.
% Code is inspired by replication codes of the paper "Booms and Banking Crises"
% by Boissay, Collard, Smets (2016, Journal of Political Economy).
% -------------------------------------------------------------------------
% INPUTS
% - ghat        [vector]   quarterly growth rate of output
% - cutoff      [scalar]   threshold for recession identification
% - minQuarters [scalar]   minimum number of quarters for a recession
% -------------------------------------------------------------------------
% OUTPUTS
% - meanDuration    [vector]   mean duration of all recessions, mild recessions, and severe recessions
% - meanOutputLoss  [vector]   mean output loss of all recessions, mild recessions, and severe recessions
% - recessionNumber [vector]   number of all recessions, severe recessions, and mild recessions
% - T               [scalar]   number of periods in ghat after removing initial recession
% =========================================================================
% Copyright (C) 2025-2026 Willi Mutschler
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
% This file is part of the replication files for the paper
% "Pruned skewed Kalman filter and smoother with application to DSGE models"
% by Gaygysyz Guljanov, Willi Mutschler, Mark Trede
% =========================================================================

%% data handling
% make sure first period is not recession
ghat = ghat(find(ghat>0,1,'first'):length(ghat));
T = length(ghat);

% building output index by simulating levels from growth rates
Output = zeros(T+1,1);
Output(1) = 1; % initialize
for t = 2:length(Output)
    Output(t) = (1+ghat(t-1))*Output(t-1);
end
Output = Output(2:end); % remove initial dummy entry to have same size as ghat

%% identify recessions
if nargin < 2 || isempty(cutoff)
    cutoff = prctile(ghat,30);
end
startRecessions = find(ghat<cutoff); % potential start of recession
recessionPeriods = zeros(T,1); % indicator of whether economy is in a recession at a given time
recessionPeriods(startRecessions) = 1; % flag the start of recessions
% extend the recession periods as long as growth is negative
for i = 1:length(startRecessions) % go through each start index
    isRecession = 1;
    tt = 1;
    while (isRecession==1) && (startRecessions(i)+tt<=T-1) % break while loop if period with nonnegative growth is encountered
        % if growth is still negative, we are still in recession. For nonnegative growth we end recession and break loop
        isRecession = (ghat(startRecessions(i)+tt)<0);
        recessionPeriods(startRecessions(i)+tt) = isRecession;
        tt = tt+1;
    end
end

%% finding recession peaks and troughs
% identify peaks (start) and troughs (end) of recessions by taking first differences on recessionPeriods
%  +1 means recessionPeriods jumped from 0 to 1 at time t+1 (i.e. recession starts next period at t+1)
%  -1 means recessionPeriods jumped from 1 to 0 at time t+1 (i.e. recession ends   next period at t+1)
recessionStartEnd = [diff(recessionPeriods);0];
recessionPeaks   = find(recessionStartEnd==+1); % peaks before recession starts
recessionTroughs = find(recessionStartEnd==-1); % troughs marking the end of recession
% match each start with its corresponding end to ensure recessionPeaks and recessionTroughs have same length (same number of recessions identified)
recessionNumber = min(length(recessionPeaks),length(recessionTroughs));
recessionPeaks   = recessionPeaks(1:recessionNumber);
recessionTroughs = recessionTroughs(1:recessionNumber);

%% discard very short recessions
% discard short recessions that last fewer than X quarters to mimick NBER's or national statistic's convention
if nargin < 3 || isempty(minQuarters)
    minQuarters = 2; % threshold how many consecutive quarters of negative growth before we call episode a recession
end
recessionPairs = [recessionPeaks((recessionTroughs-recessionPeaks)>=minQuarters) recessionTroughs((recessionTroughs-recessionPeaks)>=minQuarters)];
nRecessions = size(recessionPairs,1); % number of valid peak-trough pairs

%% compute duration and magnitude
Duration = zeros(nRecessions,1);
OutputLoss = zeros(nRecessions,1);
for i = 1:nRecessions
    % from peak ("last good quarter") to trough ("last bad quarter")
    Duration(i) = recessionPairs(i,2) - recessionPairs(i,1);
    OutputLoss(i) = 100*( Output(recessionPairs(i,2))/Output(recessionPairs(i,1)) - 1 );
end

%% grouping recessions by severity using terciles of output drop
recessionsThresholdSevere = prctile(OutputLoss,100/3);
recessionsThresholdMild   = prctile(OutputLoss,200/3);
idxSevere = OutputLoss <= recessionsThresholdSevere;
idxMild   = OutputLoss >= recessionsThresholdMild;
recessionNumber = [nRecessions  sum(idxSevere)  sum(idxMild)];
meanDuration    = [mean(Duration)   mean(Duration(idxSevere))   mean(Duration(idxMild))];
meanOutputLoss  = [mean(OutputLoss) mean(OutputLoss(idxSevere)) mean(OutputLoss(idxMild))];