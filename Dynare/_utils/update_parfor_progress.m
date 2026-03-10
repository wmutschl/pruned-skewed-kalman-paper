function update_parfor_progress(N, h, t_start)
% update_parfor_progress - GUI waitbar progress for parfor loops using DataQueue
%
% Usage:
%   q = parallel.pool.DataQueue;
%   h = waitbar(0, 'Starting...');
%   t_start = tic;
%   afterEach(q, @(~) update_parfor_progress(N, h, t_start));
%   parfor i = 1:N
%       % ... work ...
%       send(q, i);
%   end
%   delete(q);
%   close(h);
%   clear update_parfor_progress;  % reset persistent counter for next use

    persistent count
    if isempty(count)
        count = 0;
    end
    count = count + 1;
    pct = count / N;
    elapsed = toc(t_start);
    remaining = elapsed / count * (N - count);
    hrs = floor(remaining / 3600);
    mins = floor(mod(remaining, 3600) / 60);
    secs = floor(mod(remaining, 60));
    if hrs > 0
        time_str = sprintf('%dh %02dm %02ds', hrs, mins, secs);
    else
        time_str = sprintf('%02dm %02ds', mins, secs);
    end
    waitbar(pct, h, sprintf('%d/%d (%.0f%%) - Remaining: %s', count, N, pct*100, time_str));
    if count >= N
        count = [];
    end
end
