% This file is part of the replication files for the paper "Pruned Skewed
% Kalman Filter and Smoother: Pruned Skewed Kalman Filter and Smoother:
% With Applications to the Yield Curve and Asymmetric Monetary Policy Shocks"
% by Gaygysyz Guljanov, Willi Mutschler, Mark Trede

function [xmin, ...      % minimum search point of last iteration
	  fmin, ...      % function value of xmin
	  counteval, ... % number of function evaluations done
	  stopflag, ...  % stop criterion reached
	  outhist, ...   % output record with 6 columns
	  bestever ...   % struct with overall best search point
	 ] = cmaes_dsge( ...
    fitfun, ...    % name of objective/fitness function
    xstart, ...    % objective variables initial point, determines N
    sigma,  ...    % the starting value for sigma
    insigma, ...   % initial coordinate wise standard deviation(s)
    inopts, ...    % options struct, see defopts below
    varargin )     % arguments passed to objective function 
cmaVersion = '2.34'; 

% cmaes.m, Version 2.34c, last change: September, 21, 2005 
% CMAES implements an Evolution Strategy with Covariance Matrix
% Adaptation (CMA-ES) for nonlinear function minimization.  For
% introductory comments and copyright see end of file (type 'type
% cmaes').
% Adapted to DSGE models by Martin M�ller Andreasen, September 2007
%
% OPTS = CMAES returns default options. 
% OPTS = CMAES('defaults') returns default options quietly.
% OPTS = CMAES('defaults', OPTS) supplements options OPTS with default 
% options.
%
% XMIN = CMAES(FUN, X0, SIGMA, INSIGMA, [, OPTS]) locates the minimum XMIN of
% function FUN starting from column vector X0 with the initial
% coordinate wise search standard deviation SIGMA.
%
% Input arguments: 
%
%   FUN can be a function name like 'myfun' or a function handle like
%     @myfun. FUN takes as argument a column vector of size of X0 and
%     returns a scalar. An easy way to implement a hard non-linear
%     constraint is to return NaN. Then, this function evaluation is
%     not counted and a newly sampled point is tried immediately.
%
%   X0 is a column vector, or a matrix, or a string. If X0 is a matrix,
%     mean(X0, 2) is taken as initial point. If X0 is a string like
%     '2*rand(10,1)-1', the string is evaluated first.
%
% SIGMA is a scalar and determines the initial overall variance in the 
%     search. 
%
% INSIGMA is a column vector of size(X0,1), or a string
%     that can be evaluated into a column vector.  
%
%   OPTS (an optional argument) is a struct holding additional input
%     options. Valid field names and a short documentation can be
%     discovered by looking at the default options (type 'cmaes'
%     without arguments, see above). Empty or missing fields in OPTS
%     invoke the default value, i.e. OPTS needs not to have all valid
%     field names.  Capitalization does not matter and unambiguous
%     abbreviations can be used for the field names. If a string is
%     given where a numerical value is needed, the string is evaluated
%     by eval, where 'N' expands to the problem dimension
%     (==size(X0,1)). For convenience the default values can be
%     changed by editing defopts in the source code.
%
% [XMIN, FMIN, COUNTEVAL, STOPFLAG, OUTHIST, BESTEVER] = ...
%    CMAES_DSGE(FITFUN, X0, SIGMA, INSIGMA, OPTS)
% returns the best (minimal) point XMIN (found in the last
% generation); function value FMIN of XMIN; the number of needed
% function evaluations COUNTEVAL; a STOPFLAG value as cell array,
% where possible entries are 'fitness', 'tolx', 'tolupx', 'tolfun',
% 'maxfunevals', 'maxiter', 'stoptoresume' 'warnconditioncov',
% 'warnnoeffectcoord', 'warnnoeffectaxis', 'warnequalfunvals',
% 'warnequalfunvalhist', 'bug' (use e.g. any(strcmp(STOPFLAG, 'tolx'))
% or findstr(strcat(STOPFLAG, 'tolx')) for further processing); a
% history record OUTHIST with columns (1) function evaluation count,
% (2) function value, (3) axis ratio of search distribution, (4)
% maximal coordinate wise standard deviation
% (sigma*sqrt(max(diag(C)))), (5) minimal coordinate wise standard
% deviation, (6) maximal standard deviation in covariance matrix C;
% and struct BESTEVER with the overall best evaluated point x with
% function value f evaluated at evaluation count counteval.
%
% To run the code completely quietly set Display, VerboseModulo, and
% Plotting options to 0.  When OPTS.Saving==1 (default) everything is
% saved in file OPTS.SaveFileName (default 'variablescmaes.mat')
% permitting to investigate the recent result (e.g. plot with function
% plotcmaes) even while CMAES is still running (which can be quite
% useful) and to resume the search afterwards by using the resume
% option.
%
% To find the best ever evaluated point load the variables typing
% "es=load('variablescmaes')" and investigate variable es.bestever (or
% any other of the output parameters as given above). To further
% control data sampling behavior use SaveModulo option.
%
% The primary strategy parameter to play with is OPTS.PopSize, which
% can be increased from its default value.  Increasing the population
% size (by default together with the parent number OPTS.ParentNumber)
% improves global search properties in exchange to speed. Speed
% decreases, as a rule, at most linearely with increasing population
% size.
% 
% EXAMPLE
% THE OBJECTIVE FUNCTION 
% function f=frosen(x)
% f = 100*(0-x(1))^2 + 20*(0-x(2))^2+(0-x(3))^2+(0-x(4))^4;
% if  65 < x(1) & x(1) < 75
%    f = NaN;
% end
% THE CODES FOR STARTING CMAES
% x0 = 10*ones(4,1);
% insigma = 1*ones(4,1);
% sigma = 1;
% opts.TolX         = 1e-8;
% opts.TolFun       = 1e-12;
% opts.SigmaMax     = 20;
% [XMIN, FMIN, COUNTEVAL, STOPFLAG, OUTHIST, BESTEVER] = ...
%   cmaes_dsge(@frosen, x0, sigma, insigma, opts)
%
%
% ----------- Set Defaults for Input Parameters and Options -------------
% These defaults may be edited for convenience

% Input Defaults (obsolete, these are obligatory now)
definput.fitfun = 'felli'; % frosen; fcigar; see end of file for more
definput.xstart = rand(10,1); 0.50*ones(10,1);
definput.sigma = 0.3;

% Options defaults: Stopping criteria % (value of stop flag)
defopts.StopFitness  = '-Inf % stop if f(xmin) < stopfitness, minimization';
defopts.MaxFunEvals  = '3e3*((N+5)^2+popsize)  % maximal number of fevals';
defopts.MaxIter      = 'Inf  % maximal number of iterations';
defopts.StopFunEvals = 'Inf  % stop after resp. evaluation to resume later';
defopts.StopIter     = 'Inf  % stop after resp. iteration to resume later';
defopts.TolX         = '1e-12*max(insigma) % stop if x-change smaller TolX';
defopts.TolUpX       = '1e8*max(insigma) % stop if x-changes larger TolUpX';
defopts.TolFun       = '1e-12 % stop if fun-changes smaller TolFun';
defopts.StopOnWarnings = 'yes   % ''no''==''off''==0, ''on''==''yes''==1 ';
defopts.SigmaMax     = '10    % the maximum value for sigma';

% Options defaults: Other
defopts.DiffMaxChange = 'Inf  % maximal variable change(s), can be Nx1-vector';
defopts.DiffMinChange = '0    % minimal variable change(s), can be Nx1-vector';
defopts.WarnOnEqualFunctionValues = ...
    'yes   % ''no''==''off''==0, ''on''==''yes''==1 ';
defopts.LBounds = '-Inf % lower bounds, scalar or Nx1-vector'; 
defopts.UBounds = 'Inf  % upper bounds, scalar or Nx1-vector'; 

defopts.PopSize      = '(4 + floor(3*log(N)))  % population size, lambda'; 
defopts.ParentNumber = 'floor(popsize/2)     % popsize equals lambda';
defopts.RecombinationWeights = 'superlinear decrease % or linear, or equal';
defopts.Display  = 'on   % display messages like initial and final message';
defopts.Plotting = 'on   % plot while running';
defopts.VerboseModulo = '100  % >=0, messaging after every i-th iteration';
defopts.Resume  = 'no   % resume former run from SaveFile';  
defopts.Science = 'off  % off==do some additional (minor) problem capturing';
defopts.Saving  =    'on   % save data to file while running';
defopts.SaveModulo = '1    % if >1 record data less frequently after gen=100';
defopts.SaveTime   = '25   % max. percentage of time for recording data';
defopts.SaveFileName = 'variablescmaes.mat'; % file name for saving

% defopts.varopt = []; % 'for temporary and hacking purposes'; 

% ---------------------- Handling Input Parameters ----------------------

if nargin < 1 | isequal(fitfun, 'defaults') % pass default options
  if nargin < 1
    disp('Default options returned (type "help cmaes" for help).');
  end
  xmin = defopts;
  if nargin > 1 % supplement second argument with default options
    xmin = getoptions(xstart, defopts);
  end
  return;
end

input.fitfun = fitfun; % record used input
if isempty(fitfun)
  % fitfun = definput.fitfun; 
  % warning(['Objective function not determined, ''' fitfun ''' used']);
  error(['Objective function not determined']);
end
if ischar(fitfun)
  fitfun = str2func(fitfun);
end


if nargin < 2 
  xstart = [];
end

input.xstart = xstart;
if isempty(xstart)
  % xstart = definput.xstart;  % objective variables initial point
  % warning('Initial search point, and problem dimension, not determined');
  error('Initial search point, and problem dimension, not determined');
end

if nargin < 3 
  sigma = [];
end
if isa(sigma, 'struct')
  error(['Third argument SIGMA must be (or eval to) a scalar ']);
end
if isempty(sigma)
    sigma = 1;
end
input.sigma = sigma;

if nargin < 4 
  insigma = [];
end
if isa(insigma, 'struct')
  error(['Fourth argument INSIGMA must be a column vector of size(X0,1)']);
end
if isempty(insigma)
  if size(myeval(xstart),2) > 1
    insigma = std(xstart')'; 
    if any(insigma == 0)
      error(['Initial search volume is zero, choose INSIGMA or X0 appropriate']);
    end
  end
end
input.insigma = insigma;



% Compose options opts
if nargin < 5 | isempty(inopts) % no input options available
  inopts = []; 
  opts = defopts;
else
  opts = getoptions(inopts, defopts);
end

% ------------------------ Initialization -------------------------------

% Handle resuming of old run
flgresume = myevalbool(opts.Resume);
if ~flgresume % not resuming a former run
  % Assign settings from input parameters and options for myeval...
  xmean = mean(myeval(xstart), 2); % in case of xstart is a population 
  N = size(xmean, 1); numberofvariables = N; 
  popsize = myeval(opts.PopSize); lambda = popsize;
  insigma = myeval(insigma);
  if all(size(insigma) == [N 2]) 
    insigma = 0.5 * (insigma(:,2) - insigma(:,1));
  end
else % flgresume is true, do resume former run
  tmp = whos('-file', opts.SaveFileName);
  for i = 1:length(tmp)
    if strcmp(tmp(i).name, 'localopts');
      error('Saved variables include variable "localopts", please remove');
    end
  end
  localopts = opts; % keep stopping and display options
  load(opts.SaveFileName); 
  flgresume = 1;
  
  % Overwrite old stopping and display options
  opts.StopFitness = localopts.StopFitness; 
  opts.MaxFunEvals = localopts.MaxFunEvals;
  opts.MaxIter = localopts.MaxIter; 
  opts.StopFunEvals = localopts.StopFunEvals; 
  opts.StopIter = localopts.StopIter;  
  opts.TolX = localopts.TolX;
  opts.TolUpX = localopts.TolUpX;
  opts.TolFun = localopts.TolFun;
  opts.StopOnWarnings = localopts.StopOnWarnings; 
  opts.SigmaMax = localopts.SigmaMax;
  opts.Display = localopts.Display;
  opts.Plotting = localopts.Plotting;
  opts.VerboseModulo = localopts.VerboseModulo;
  opts.Saving = localopts.Saving;
  opts.SaveModulo = localopts.SaveModulo;
  opts.SaveTime = localopts.SaveTime;
  clear localopts; % otherwise localopts would be overwritten during load
end
  
% Evaluate options
stopFitness = myeval(opts.StopFitness); 
stopMaxFunEvals = myeval(opts.MaxFunEvals);  
stopMaxIter = myeval(opts.MaxIter);  
stopFunEvals = myeval(opts.StopFunEvals);  
stopIter = myeval(opts.StopIter);  
stopTolX = myeval(opts.TolX);
stopTolUpX = myeval(opts.TolUpX);
stopTolFun = myeval(opts.TolFun);
stopOnWarnings = myevalbool(opts.StopOnWarnings); 
SigmaMax = myeval(opts.SigmaMax);
flgWarnOnEqualFunctionValues = myevalbool(opts.WarnOnEqualFunctionValues);
flgdisplay = myevalbool(opts.Display);
flgplotting = myevalbool(opts.Plotting);
verbosemodulo = myeval(opts.VerboseModulo);
flgscience = myevalbool(opts.Science);
flgsaving = myevalbool(opts.Saving);
savemodulo = myeval(opts.SaveModulo);
savetime = myeval(opts.SaveTime);

if (isfinite(stopFunEvals) | isfinite(stopIter)) & ~flgsaving
  warning('To resume later the saving option needs to be set');
end

% Do more checking and initialization 
if ~flgresume
  maxdx = myeval(opts.DiffMaxChange); % maximal sensible variable change
  mindx = myeval(opts.DiffMinChange); % minimal sensible variable change 
				      % can both also be defined as Nx1 vectors
  lbounds = myeval(opts.LBounds);		     
  ubounds = myeval(opts.UBounds);
  if isempty(insigma) % last chance to set insigma
    if all(lbounds > -Inf) & all(ubounds < Inf)
      if any(lbounds>=ubounds)
	error('upper bound must be greater than lower bound');
      end
      insigma = 0.3*(ubounds-lbounds);
      stopTolX = myeval(opts.TolX);  % reevaluate these
      stopTolUpX = myeval(opts.TolUpX);
    else
      error(['Initial step sizes (INSIGMA) not determined']);
    end
  end

  % Check all vector sizes
  if size(xmean, 2) > 1 | size(xmean,1) ~= N
    error(['intial search point should be a column vector of size ' ...
	   num2str(N)]);
  elseif ~(all(size(insigma) == [1 1]) | all(size(insigma) == [N 1]))
    error(['input parameter INSIGMA should be (or eval to) a scalar '...
	   'or a column vector of size ' num2str(N)] );
  elseif size(stopTolX, 2) > 1 | ~ismember(size(stopTolX, 1), [1 N])
    error(['option TolX should be (or eval to) a scalar '...
	   'or a column vector of size ' num2str(N)] );
  elseif size(stopTolUpX, 2) > 1 | ~ismember(size(stopTolUpX, 1), [1 N])
    error(['option TolUpX should be (or eval to) a scalar '...
	   'or a column vector of size ' num2str(N)] );
  elseif size(maxdx, 2) > 1 | ~ismember(size(maxdx, 1), [1 N])
    error(['option DiffMaxChange should be (or eval to) a scalar '...
	   'or a column vector of size ' num2str(N)] );
  elseif size(mindx, 2) > 1 | ~ismember(size(mindx, 1), [1 N])
    error(['option DiffMinChange should be (or eval to) a scalar '...
	   'or a column vector of size ' num2str(N)] );
  elseif size(lbounds, 2) > 1 | ~ismember(size(lbounds, 1), [1 N])
    error(['option lbounds should be (or eval to) a scalar '...
	   'or a column vector of size ' num2str(N)] );
  elseif size(ubounds, 2) > 1 | ~ismember(size(ubounds, 1), [1 N])
    error(['option ubounds should be (or eval to) a scalar '...
	   'or a column vector of size ' num2str(N)] );
  end
  
  % Strategy internal parameter setting: Selection  
  lambda = myeval(opts.PopSize);  % population size, offspring number
  mu = myeval(opts.ParentNumber); % number of parents/points for recombination
  if strncmp(lower(opts.RecombinationWeights), 'equal', 3)
    weights = ones(mu,1); % (mu_I,lambda)-CMA-ES
  elseif strncmp(lower(opts.RecombinationWeights), 'linear', 3)
    weights = mu+1-(1:mu)';
  elseif strncmp(lower(opts.RecombinationWeights), 'superlinear', 3)
    weights = log(mu+1)-log(1:mu)'; % muXone array for weighted recombination
  else
    error(['Recombination weights to be "' opts.RecombinationWeights ...
	   '" is not implemented']);
  end
  mueff=sum(weights)^2/sum(weights.^2); % variance-effective size of mu
  if mueff == lambda
    error(['Combination of values for PopSize, ParentNumber and ' ...
	  ' and RecombinationWeights is not reasonable']);
  end
  
  % Strategy internal parameter setting: Adaptation
  cc = 4/(N+4);         % time constant for cumulation for covariance matrix
  cs = (mueff+2)/(N+mueff+3); % t-const for cumulation for step size control
  mucov = mueff;   % size of mu used for calculating learning rate ccov
  ccov = (1/mucov) * 2/(N+1.41)^2 ... % learning rate for covariance matrix
	 + (1-1/mucov) * min(1,(2*mueff-1)/((N+2)^2+mueff)); 
  % ||ps|| is close to sqrt(mueff/N) for mueff large on linear fitness
  damps = ... % damping for step size control, usually close to one 
      (1 + 2*max(0,sqrt((mueff-1)/(N+1))-1)) ... % limit sigma increase
      * max(0.3, ... % reduce damps, if max. iteration number is small
	  1 - N/min(stopMaxIter,stopMaxFunEvals/lambda)) + cs; 

  %qqq hacking of a different parameter setting, e.g. for ccov or damps,
  % can be done here. 

  % Initialize dynamic internal strategy parameters
  if any(insigma <= 0) 
    error(['Initial search volume (INSIGMA) must have elements greater than zero']);
  end
  if max(insigma)/min(insigma) > 1e6
    error(['Initial search volume (INSIGMA) badly conditioned']);
  end
  pc = zeros(N,1); ps = zeros(N,1);  % evolution paths for C and sigma
  if length(insigma) == 1
    insigma = insigma * ones(N,1) ;
  end
  B = eye(N);                        % B defines the coordinate system
  D = diag(insigma);                 % diagonal matrix D defines the scaling
  BD = B*D;                          % for speed up only
  C = BD*(BD)';                      % covariance matrix
  fitness.hist=NaN*ones(1,10+ceil(3*10*N/lambda)); % history of fitness values
  fitness.histsel=NaN*ones(1,10+ceil(3*10*N/lambda)); % history of fitness values

  % Initialize boundary handling
  bnd.isactive = any(lbounds > -Inf) | any(ubounds < Inf); 
  if bnd.isactive
    if any(lbounds>ubounds)
      error('lower bound found to be greater than upper bound');
    end
    [xmean ti] = xintobounds(xmean, lbounds, ubounds); % just in case
    if any(ti)
      warning('Initial point was out of bounds, corrected');
    end
    bnd.weights = zeros(N,1);         % weights for bound penalty
    bnd.scale = diag(C)/mean(diag(C));
    bnd.isbounded = (lbounds > -Inf) | (ubounds < Inf);
    if length(bnd.isbounded) == 1
      bnd.isbounded = bnd.isbounded * ones(N,1);
    end
    maxdx = min(maxdx, (ubounds - lbounds)/2);
    if any(sigma*sqrt(diag(C)) > maxdx)
      fac = min(maxdx ./ sqrt(diag(C)))/sigma;
      insigma = min(maxdx ./ sqrt(diag(C)));
      warning(['Initial INSIGMA multiplied by the factor ' num2str(fac) ...
	       ', because it was larger than half' ...
	       ' of one of the boundary intervals']);
    end
    idx = (lbounds > -Inf) & (ubounds < Inf);
    dd = diag(C);
    if any(5*sigma*sqrt(dd(idx)) < ubounds(idx) - lbounds(idx))
      warning(['Initial INSIGMA is, in at least one coordinate, ' ...
	       'much smaller than the '...
	       'given boundary intervals. For reasonable ' ...
	       'global search performance INSIGMA should be ' ...
	       'between 0.2 and 0.5 of the bounded interval in ' ...
	       'each coordinate. If all coordinates have ' ... 
	       'lower and upper bounds INSIGMA can be empty']);
    end
    bnd.dfithist = 1;              % delta fit for setting weights
    bnd.aridxpoints = [];          % remember complete outside points
    bnd.arfitness = [];            % and their fitness
    bnd.validfitval = 0;
    bnd.iniphase = 1;
  end

  % ooo initial feval, for output only, (un-)comment if required
  counteval = 0; 
  fitness.hist(1)=feval(fitfun, xmean, varargin{:}); 
  fitness.histsel(1)=fitness.hist(1);
  counteval = counteval + 1;
                                         
  % Initialize further constants
  randn('state', sum(100*clock));     % random number generator state
  startseed = randn('state');         % for saving purpose
  chiN=N^0.5*(1-1/(4*N)+1/(21*N^2));  % expectation of 
				      %   ||N(0,I)|| == norm(randn(N,1))
  weights = weights/sum(weights);     % normalize recombination weights array
  
  % Initialize records and output
  time.t0 = clock;
  outhist = [0 fitness.hist(1) max(diag(D))/min(diag(D)) ...
	     sigma*sqrt(max(diag(C))) ...
	     sigma*sqrt(min(diag(C))) sqrt(max(diag(C)))];
  outt0=clock; outetime=0; 
  out.x = 0;                   
  out.y1=[fitness.hist(1) sigma max(diag(D))/min(diag(D)) ...
	  sigma*[max(diag(D)) min(diag(D))]];
  out.y2=xmean'; out.y2a=xmean';
  out.y3=sigma*sqrt(diag(C))';
  out.y4=sort(diag(D))'; 
  outiter = 0;

  countiter = 0;

else % resume is on
  if flgdisplay
    disp(['  resumed from ' opts.SaveFileName ]); 
  end
  if counteval >= stopMaxFunEvals 
    error(['MaxFunEvals exceeded, use StopFunEvals as stopping ' ...
	  'criterion before resume']);
  end
  if countiter >= stopMaxIter 
    error(['MaxIter exceeded, use StopIter as stopping criterion ' ...
	  'before resume']);
  end
  
end % else, if ~flgresume

% Display initial message
if flgdisplay
  if mu < 8
    disp(['  n=' num2str(N) ': (' num2str(mu) ',' ...
	    num2str(lambda) ')-CMA-ES(w=[' ...
	    num2str(100*weights(1:end-1)','%.0f ') ... 
	    	    num2str(100*weights(end)','%.0f') ']%, ' ...
	    'mu_eff=' num2str(mueff,'%.1f') ...
	    ') on function ' ...
	    func2str(fitfun) ]);
  else
    disp(['  n=' num2str(N) ': (' num2str(mu) ',' ...
	    num2str(lambda) ')-CMA-ES (w=[' ...
	    num2str(100*weights(1:2)','%.2g ') ...
	    num2str(100*weights(3)','%.2g') '...' ...
	    num2str(100*weights(end-1:end)',' %.2g') ']%, ' ...
	    'mu_eff=' num2str(mueff,'%.1f') ...
	    ') on function ' ...
	    func2str(fitfun)]);
  end
end

% -------------------- Generation Loop --------------------------------
stopflag = {};
while isempty(stopflag)
  countiter = countiter + 1; 
  
  % Generate and evaluate lambda offspring
 
  for k=1:lambda,
    fitness.raw(k) = NaN; 
    tries = 0;
    % Resample, until fitness is not NaN
    while isnan(fitness.raw(k))
      arz(:,k) = randn(N,1);
      arx(:,k) = xmean + sigma * (BD * arz(:,k));                % Eq. (1)

      % You may handle constraints here. You may either resample
      % arz(:,k) and/or multiply it with a factor between -1 and 1
      % (the latter will decrease the overall step size) and
      % recalculate arx accordingly. Do not change arx or arz in any
      % other way.
 
      if ~bnd.isactive
        arxvalid(:,k) = arx(:,k);
      else
        arxvalid(:,k) = xintobounds(arx(:,k), lbounds, ubounds);
      end
      % You may handle constraints here.  You may copy and alter
      % (columns of) arxvalid(:,k) only for the evaluation of the
      % fitness function. arx and arxvalid should not be changed.
      fitness.raw(k) = feval(fitfun, arxvalid(:,k), varargin{:}); 
      tries = tries + 1;
      if tries==500*lambda 
      	warning(['Many NaN objective function values at evaluation ' ...
		num2str(counteval)]);
        sigma = 0.9*sigma;  
        tries = 0;
        disp(['sigma reduced by 10%. Sigma = ', num2str(sigma)]);
      end
    end
    % NaN are not counted
    counteval = counteval + 1;
  end

  fitness.sel = fitness.raw; 

  % ----- handle boundaries -----
  if 1 < 3 & bnd.isactive
    % Get delta fitness values
    val = myprctile(fitness.raw, [25 75]);
    val = (val(2) - val(1)) / N / mean(diag(C)) / sigma^2;
    %val = (myprctile(fitness.raw, 75) - myprctile(fitness.raw, 25)) ...
    %    / N / mean(diag(C)) / sigma^2;
    % Catch non-sensible values 
    if ~isfinite(val)
      %warning('Non-finite fitness range');
      val = max(bnd.dfithist);  
    elseif val == 0 % happens if all points are out of bounds
      val = min(bnd.dfithist(bnd.dfithist>0)); 
    elseif bnd.validfitval == 0 % first sensible val
      bnd.dfithist = [];
      bnd.validfitval = 1;
    end

    % Keep delta fitness values
    if length(bnd.dfithist) < 20+(3*N)/lambda
      bnd.dfithist = [bnd.dfithist val];
    else
      bnd.dfithist = [bnd.dfithist(2:end) val];
    end

    % Scale weights anew, bias scaling to unity
    if 1 < 3
      bnd.weights = bnd.scale .* bnd.weights;  % reset scaling
      bnd.scale = exp(0.1*mean(log(diag(C)))) * diag(C).^0.9;
      bnd.scale = bnd.scale / exp(mean(log(bnd.scale))); % prod is 1 initially
      bnd.weights = bnd.weights ./ bnd.scale; 
    end
    [tx ti]  = xintobounds(xmean, lbounds, ubounds);

    % Set initial weights
    if bnd.iniphase 
      if any(ti) 
        bnd.weights(find(bnd.isbounded)) = ...
          2.0002 * median(bnd.dfithist) ./ bnd.scale(find(bnd.isbounded));
	if bnd.validfitval & countiter > 2
          bnd.iniphase = 0;
        end
      end
    end

    % Increase/decrease weights
    if  1 < 3 & any(ti) % any coordinate of xmean out of bounds
      % judge distance of xmean to boundary
      tx = xmean - tx;
      idx = (ti ~= 0 & abs(tx) > 3*max(1,sqrt(N)/mueff) ... 
	     * sigma*sqrt(diag(C))) ;
      if ~isempty(idx) % increase
	bnd.weights(idx) = 1.1^(max(1, mueff/10/N)) * bnd.weights(idx); 
      end
    end

    % Assigned penalized fitness
    bnd.arpenalty = bnd.weights' * (arxvalid - arx).^2; 
    fitness.sel = fitness.raw + bnd.arpenalty;

  end % handle boundaries
  % ----- end handle boundaries -----
  
  % Sort by fitness 
  [fitness.raw, fitness.idx] = sort(fitness.raw);  
  [fitness.sel, fitness.idxsel] = sort(fitness.sel);   % minimization
  fitness.hist(2:end) = fitness.hist(1:end-1);    % record short history of
  fitness.hist(1) = fitness.raw(1);               % best fitness values
  fitness.histsel(2:end) = fitness.histsel(1:end-1);    % record short history of
  fitness.histsel(1) = fitness.sel(1);               % best fitness values

  % Calculate new xmean, this is selection and recombination 
  xold = xmean; % for speed up of Eq. (2) and (3)
  xmean = arx(:,fitness.idxsel(1:mu))*weights; 
  zmean = arz(:,fitness.idxsel(1:mu))*weights;%==D^-1*B'*(xmean-xold)/sigma
  % fmean = feval(fitfun, xintobounds(xmean, lbounds, ubounds), varargin{:});
  % counteval = counteval + 1;
  
  % Cumulation: update evolution paths
  ps = (1-cs)*ps + (sqrt(cs*(2-cs)*mueff)) * (B*zmean);          % Eq. (4)
  hsig = norm(ps)/sqrt(1-(1-cs)^(2*countiter))/chiN < 1.5 + 1/(N-0.5);
  pc = (1-cc)*pc ...
        + hsig*(sqrt(cc*(2-cc)*mueff)/sigma) * (xmean-xold);     % Eq. (2)
  if hsig == 0
    %disp([num2str(countiter) ' ' num2str(counteval) ' pc update stalled']);
  end

  % Adapt covariance matrix
  if ccov > 0                                                    % Eq. (3)
    C = (1-ccov+(1-hsig)*ccov*cc*(2-cc)/mucov) * C ... % regard old matrix 
        + ccov * (1/mucov) * pc*pc' ...                % plus rank one update
        + ccov * (1-1/mucov) ...                       % plus rank mu update
          * sigma^-2 * (arx(:,fitness.idxsel(1:mu))-repmat(xold,1,mu)) ...
          * diag(weights) * (arx(:,fitness.idxsel(1:mu))-repmat(xold,1,mu))';
  end
  
  if 1 < 2 & ~flgscience 
    % remove momentum in ps, if ps is large and fitness is getting worse.
    % this should rarely happen. 
    % this might be questionable in dynamic environments
    if sum(ps.^2)/N > 1.5 + 10*(2/N)^.5 & ...
        fitness.histsel(1) > max(fitness.histsel(2:3))
      ps = ps * sqrt(N*(1+max(0,log(sum(ps.^2)/N))) / sum(ps.^2));
      if flgdisplay
        disp(['Momentum in ps removed at [niter neval]=' ...
              num2str([countiter counteval]) ']']);
      end
    end
  end

  % Adapt sigma 
  sigma = min(sigma * exp((norm(ps)/chiN - 1)*cs/damps),SigmaMax);             % Eq. (5)

  % Update B and D from C
  if ccov > 0 & mod(countiter, 1/ccov/N/10) < 1
    C=triu(C)+triu(C,1)'; % enforce symmetry
    [B,D] = eig(C);       % eigen decomposition, B==normalized eigenvectors
    % limit condition of C to 1e14 + 1
    if min(diag(D)) <= 0
	if stopOnWarnings
	  stopflag(end+1) = {'warnconditioncov'};
	else
	  warning(['Iteration ' num2str(countiter) ...
		   ': Eigenvalue (smaller) zero']);
	  D(D<0) = 0;
	  tmp = max(diag(D))/1e14;
	  C = C + tmp*eye(N); D = D + tmp*eye(N); 
	end
    end
    if max(diag(D)) > 1e14*min(diag(D)) 
	if stopOnWarnings
	  stopflag(end+1) = {'warnconditioncov'};
	else
	  warning(['Iteration ' num2str(countiter) ': condition of C ' ...
		   'at upper limit' ]);
	  tmp = max(diag(D))/1e14 - min(diag(D));
	  C = C + tmp*eye(N); D = D + tmp*eye(N); 
	end
    end
    D = diag(sqrt(diag(D))); % D contains standard deviations now
    % D = D / prod(diag(D))^(1/N);  C = C / prod(diag(D))^(2/N);
    BD = B*D; % for speed up only
  end % if mod

  % ----- numerical error management -----
  % Adjust maximal coordinate axis deviations
  if any(sigma*sqrt(diag(C)) > maxdx)
    sigma = min(maxdx ./ sqrt(diag(C)));
    %warning(['Iteration ' num2str(countiter) ': coordinate axis std ' ...
    %         'deviation at upper limit of ' num2str(maxdx)]);
    % stopflag(end+1) = {'maxcoorddev'};
  end
  % Adjust minimal coordinate axis deviations
  while any(sigma*sqrt(diag(C)) < mindx)
    sigma = max(mindx ./ sqrt(diag(C))) * exp(0.05+cs/damps); 
    %warning(['Iteration ' num2str(countiter) ': coordinate axis std ' ...
    %         'deviation at lower limit of ' num2str(mindx)]);
    % stopflag(end+1) = {'mincoorddev'};;
  end
  % Adjust too low coordinate axis deviations
  if any(xmean == xmean + 0.2*sigma*sqrt(diag(C))) 
    if stopOnWarnings
	stopflag(end+1) = {'warnnoeffectcoord'};
    else
      warning(['Iteration ' num2str(countiter) ': coordinate axis std ' ...
	       'deviation too low' ]);
	C = C + ccov * diag(diag(C) .* ...
			    (xmean == xmean + 0.2*sigma*sqrt(diag(C))));
	sigma = sigma * exp(0.05+cs/damps); 
    end
  end
  % Adjust step size in case of (numerical) precision problem 
  if all(xmean == xmean ...
	    + 0.1*sigma*BD(:,1+floor(mod(countiter,N))))
    i = 1+floor(mod(countiter,N));
    if stopOnWarnings
	stopflag(end+1) = {'warnnoeffectaxis'};
    else
      warning(['Iteration ' num2str(countiter) ...
	       ': main axis standard deviation ' ...
	       num2str(sigma*D(i,i)) ' has no effect' ]);
	sigma = sigma * exp(0.2+cs/damps); 
    end
  end
  % Adjust step size in case of equal function values (flat fitness)
  if fitness.sel(1) == fitness.sel(1+ceil(0.1+lambda/4))
    if flgWarnOnEqualFunctionValues & stopOnWarnings
	stopflag(end+1) = {'warnequalfunvals'};
    else
      if flgWarnOnEqualFunctionValues
	warning(['Iteration ' num2str(countiter) ...
		 ': equal function values f=' num2str(fitness.sel(1)) ...
		 ' at maximal main axis sigma ' ...
		 num2str(sigma*max(diag(D)))]);
      end
      sigma = sigma * exp(0.2+cs/damps); 
    end
  end
  % Adjust step size in case of equal function values
  if countiter > 2 & myrange([fitness.hist fitness.sel(1)]) == 0  
    if stopOnWarnings
	stopflag(end+1) = {'warnequalfunvalhist'};
    else
      warning(['Iteration ' num2str(countiter) ...
	       ': equal function values in history at maximal main ' ...
	       'axis sigma ' num2str(sigma*max(diag(D)))]);
	sigma = sigma * exp(0.2+cs/damps); 
    end
  end
  % Align scales of sigma and C for nicer output
  if 11 < 2 & (sigma > 1e10*max(diag(D)) | 1e10*sigma < min(diag(D)))
    fac = sqrt(sigma/median(diag(D)));
    sigma = sigma/fac;
    pc = fac * pc;
    C = fac^2 * C;
    D = fac * D;
    BD = fac * BD;
  end
    
  % ----- end numerical error management -----
  
  % Keep overall best solution
  if countiter == 1
    if fitness.hist(2) < fitness.hist(1)
	bestever.x = xold;
	bestever.f = fitness.hist(2);
	bestever.counteval = 1;
    else
	bestever.x = arxvalid(:, fitness.idx(1));
	bestever.f = fitness.hist(1);
	bestever.counteval = counteval + fitness.idx(1) - lambda;
    end
  elseif fitness.hist(1) < bestever.f
    bestever.x = arxvalid(:, fitness.idx(1));
    bestever.f = fitness.hist(1);
    bestever.counteval = counteval + fitness.idx(1) - lambda;
  end

  % Set stop flag
  if fitness.raw(1) <= stopFitness stopflag(end+1) = {'fitness'}; end
  if counteval >= stopMaxFunEvals stopflag(end+1) = {'maxfunevals'}; end
  if countiter >= stopMaxIter stopflag(end+1) = {'maxiter'}; end
  if all(sigma*(max(abs(pc), sqrt(diag(C)))) < stopTolX) 
    stopflag(end+1) = {'tolx'};
  end
  if any(sigma*sqrt(diag(C)) > stopTolUpX) 
    stopflag(end+1) = {'tolupx'};
  end
  if sigma*max(D) == 0  % should never happen
    stopflag(end+1) = {'bug'};
  end
  if countiter > 2 & myrange([fitness.sel fitness.hist]) < stopTolFun 
    stopflag(end+1) = {'tolfun'};
  end
  if counteval >= stopFunEvals | countiter >= stopIter
    stopflag(end+1) = {'stoptoresume'};
    if length(stopflag) == 1 & flgsaving == 0
      error('To resume later the saving option needs to be set');
    end
  end

  % ----- output generation -----
  if verbosemodulo > 0 & isfinite(verbosemodulo)
    if countiter == 1 | mod(countiter, 10*verbosemodulo) < 1 
      disp(['Iterat, #Fevals:   Function Value    (median,worst) ' ...
	    '|Axis Ratio|' ...
	    'idx:Min SD idx:Max SD']); 
    end
    if mod(countiter, verbosemodulo) < 1 ...
	  | (countiter < 3 & verbosemodulo > 0 & isfinite(verbosemodulo))
      [minstd minstdidx] = min(sigma*sqrt(diag(C)));
      [maxstd maxstdidx] = max(sigma*sqrt(diag(C)));
      % format display nicely
      disp([repmat(' ',1,4-floor(log10(countiter))) ...
	    num2str(countiter) ' , ' ...
	    repmat(' ',1,5-floor(log10(counteval))) ...
	    num2str(counteval) ' : ' ...
            num2str(fitness.hist(1), '%.13e') ...
	    ' +(' num2str(median(fitness.raw)-fitness.hist(1), '%.0e ') ...
	    ',' num2str(max(fitness.raw)-fitness.hist(1), '%.0e ') ...
	    ') | ' ...
	    num2str(max(diag(D))/min(diag(D)), '%4.2e') ' | ' ...
	    repmat(' ',1,1-floor(log10(minstdidx))) num2str(minstdidx) ':' ...
	    num2str(minstd, ' %.1e') ' ' ...
	    repmat(' ',1,1-floor(log10(maxstdidx))) num2str(maxstdidx) ':' ...
	    num2str(maxstd, ' %.1e')]);
    end
  end

  % measure time for recording data
  if countiter < 3
    time.c = 0.5;
    time.nonoutput = 0;
    time.recording = 0;
    time.saving  = 0.5; % first saving after 10 seconds
    time.plotting = 0;
  else
    time.c = min(1, time.nonoutput/3 + 1e-9); % set backward horizon
    time.c = max(1e-5, 1/countiter); % mean over all or 1e-5
  end
  % get average time per iteration
  time.t1 = clock;
  time.act = max(0,etime(time.t1, time.t0));
  time.nonoutput = (1-time.c) * time.nonoutput ...
      + time.c * time.act; 

  time.recording = (1-time.c) * time.recording;
  time.saving  = (1-time.c) * time.saving;
  time.plotting = (1-time.c) * time.plotting;
  
  % record output data, concerning time issues
  if countiter < 1e2 | ~isempty(stopflag) | ...
	countiter >= outiter + savemodulo
    outiter = countiter; 
      % Compose output argument No 5
      outhist = [outhist; [counteval fitness.hist(1) ...
		    max(diag(D))/min(diag(D)) ...
		    sigma*sqrt(max(diag(C))) sigma*sqrt(min(diag(C))) ...
		    sqrt(max(diag(C)))]];

    if (flgsaving | flgplotting)
      out.x = [out.x counteval];
      out.y1 = [out.y1; [fitness.raw(1) sigma max(diag(D))/min(diag(D))] ...
                sigma*[max(diag(D)) min(diag(D))]];
      out.y2 = [out.y2; xmean'];
      out.y2a = [out.y2a; (arx(:,fitness.idx(1)))'];
      %out.y2 = [out.y2; (arxvalid(:,fitness.idx(1)))'];
      %out.y2a = [out.y2a; xmean'];
      out.y3 = [out.y3; sigma*sqrt(diag(C))'];
      out.y4 = [out.y4; sort(diag(D))']; 
    end % end flgsaving or flgplotting
    
    % get average time for recording data
    time.t2 = clock;
    time.recording = time.recording + time.c * max(0,etime(time.t2, time.t1)); 
    
    if flgplotting & countiter > 1
      if ~isempty(stopflag) || ...
	  time.plotting < 0.2 * time.nonoutput
	outplot(out); % outplot defined below
	if countiter > 3
	  time.plotting = time.plotting + time.c * max(0,etime(clock, time.t2)); 
	end
      end
    end
    if countiter > 100 & ...
	  time.recording > savetime * (time.nonoutput+time.recording) / 100
      savemodulo = savemodulo + 1;
      % disp('++savemodulo'); %qqq
    end
  end % if output

  % save everything
  time.t3 = clock;
  if ~isempty(stopflag) | time.saving < 0.05 * time.nonoutput 
    xmin = arxvalid(:, fitness.idx(1));
    fmin = fitness.raw(1);
    if flgsaving & countiter > 2
      save(opts.SaveFileName); % for inspection and possible restart	
      time.saving = time.saving + time.c * max(0,etime(clock, time.t3)); 
    end
  end
  time.t0 = clock;

  % ----- end output generation -----
  
end % while, end generation loop

% -------------------- Final Procedures -------------------------------

% Evaluate xmean and return best recent point in xmin
fmin = fitness.raw(1);
xmin = arxvalid(:, fitness.idx(1)); % Return best point of last generation.
if length(stopflag) > sum(strcmp(stopflag, 'stoptoresume')) % final stopping
  fmean = feval(fitfun, xintobounds(xmean, lbounds, ubounds), varargin{:});
  counteval = counteval + 1;
  if fmean < fitness.raw(1)
    fmin = fmean;
    xmin = xintobounds(xmean, lbounds, ubounds); % Return xmean as best point
  end
end

% Save everything and display final message
if flgsaving
  save(opts.SaveFileName);    % for inspection and possible restart
  message = [' (saved to ' opts.SaveFileName ')'];
else
  message = [];
end

if flgdisplay
  disp(['#Fevals:   f(returned x)   |    bestever.f     | stopflag' ...
        message]);
  strstop = strcat(stopflag, '.');
  disp([repmat(' ',1,6-floor(log10(counteval))) ...
        num2str(counteval, '%6.0f') ': ' num2str(fmin, '%.11e') ' | ' ...
        num2str(bestever.f, '%.11e') ' | ' ...
	strstop{1:end}]);
  if exist('sfile', 'var') disp(['Results saved in ' sfile]); 
  end
end

% ---------------------------------------------------------------  
% ---------------------------------------------------------------  
function [x, idx] = xintobounds(x, lbounds, ubounds)
  if ~isempty(lbounds)
    idx = x < lbounds;
    if length(lbounds) == 1
      x(idx) = lbounds;
    else
      x(idx) = lbounds(idx);
    end
  end
  if ~isempty(ubounds)
    idx2 = x > ubounds;
    if length(ubounds) == 1
      x(idx2) = ubounds;
    else
      x(idx2) = ubounds(idx2);
    end
  end
  idx = idx2-idx;
  
% ---------------------------------------------------------------  
% ---------------------------------------------------------------  
function opts=getoptions(inopts, defopts)
% OPTS = GETOPTIONS(INOPTS, DEFOPTS) handles an arbitrary number of
% optional arguments to a function. The given arguments are collected
% in the struct INOPTS.  GETOPTIONS matches INOPTS with a default
% options struct DEFOPTS and returns the merge OPTS.  Empty or missing
% fields in INOPTS invoke the default value.  Fieldnames in INOPTS can
% be abbreviated.
%
% The returned struct OPTS is first assigned to DEFOPTS. Then any
% field value in OPTS is replaced by the respective field value of
% INOPTS if (1) the field unambiguously (case-insensitive) matches
% with the fieldname in INOPTS (cut down to the length of the INOPTS
% fieldname) and (2) the field is not empty.
%
% Example:
%   In the source-code of the function that needs optional
%   arguments, the last argument is the struct of optional
%   arguments:
%
%   function results = myfunction(mandatory_arg, inopts)
%     % Define four default options
%     defopts.PopulationSize = 200;
%     defopts.ParentNumber = 50;
%     defopts.MaxIterations = 1e6;
%     defopts.MaxSigma = 1;
%  
%     % merge default options with input options
%     opts = getoptions(inopts, defopts);
%
%     % Thats it! From now on the values in opts can be used
%     for i = 1:opts.PopulationSize
%       % do whatever
%       if sigma > opts.MaxSigma
%         % do whatever
%       end
%     end
%   
%   For calling the function myfunction with default options:
%   myfunction(argument1, []);
%   For calling the function myfunction with modified options:
%   opt.pop = 100; % redefine PopulationSize option
%   opt.PAR = 10;  % redefine ParentNumber option
%   opt.maxiter = 2; % opt.max=2 is ambiguous and would result in an error 
%   myfunction(argument1, opt);

%
% 04/07/19: Entries can be structs itself leading to a recursive
%           call to getoptions. 
%

if nargin < 2 | isempty(defopts) % no default options available
  opts=inopts;
  return;
elseif isempty(inopts) % empty inopts invoke default options
  opts = defopts;
  return;
elseif ~isstruct(defopts) % handle a single option value
  if isempty(inopts) 
    opts = defopts;
  elseif ~isstruct(inopts)
    opts = inopts;
  else
    error('Input options are a struct, while default options are not');
  end
  return;
elseif ~isstruct(inopts) % no valid input options
  error('The options need to be a struct or empty');
end

  opts = defopts; % start from defopts 
  % if necessary overwrite opts fields by inopts values
  defnames = fieldnames(defopts);
  idxmatched = []; % indices of defopts that already matched
  for name = fieldnames(inopts)'
    name = name{1}; % name of i-th inopts-field
    idx = strncmp(lower(defnames), lower(name), length(name));
    if sum(idx) > 1
      error(['option "' name '" is not an unambigous abbreviation. ' ...
	     'Use opts=RMFIELD(opts, ''' name, ...
	     ''') to remove the field from the struct.']);
    end
    if sum(idx) == 1
      defname  = defnames{find(idx)}; 
      if ismember(find(idx), idxmatched)
	error(['input options match more than ones with "' ...
	       defname '". ' ...
	       'Use opts=RMFIELD(opts, ''' name, ...
	       ''') to remove the field from the struct.']);
      end
      idxmatched = [idxmatched find(idx)];
      val = getfield(inopts, name);
      % next line can replace previous line from MATLAB version 6.5.0 on
      % val = inopts.(name);
      if isstruct(val) % valid syntax only from version 6.5.0
	opts = setfield(opts, defname, ...
	    getoptions(val, getfield(defopts, defname))); 
      elseif isstruct(getfield(defopts, defname)) 
      % next three lines can replace previous three lines from MATLAB 
      % version 6.5.0 on
      %   opts.(defname) = ...
      %      getoptions(val, defopts.(defname)); 
      % elseif isstruct(defopts.(defname)) 
	warning(['option "' name '" disregarded (must be struct)']); 
      elseif ~isempty(val) % empty value: do nothing, i.e. stick to default
	opts = setfield(opts, defnames{find(idx)}, val);
	% next line can replace previous line from MATLAB version 6.5.0 on
	% opts.(defname) = inopts.(name); 
      end
    else
      warning(['option "' name '" disregarded (unknown field name)']);
    end
  end

% ---------------------------------------------------------------  
% ---------------------------------------------------------------  
function res=myeval(s)
  if ischar(s)
    res = evalin('caller', s);
  else
    res = s;
  end
  
% ---------------------------------------------------------------  
% ---------------------------------------------------------------  
function res=myevalbool(s)
  if ~ischar(s) % s may not and cannot be empty
    res = s;
  else % evaluation string s
    if strncmp(lower(s), 'yes', 3) | strncmp(lower(s), 'on', 2) ...
	  | strncmp(lower(s), 'true', 4) | strncmp(s, '1 ', 2)
      res = 1;
    elseif strncmp(lower(s), 'no', 2) | strncmp(lower(s), 'off', 3) ...
	  | strncmp(lower(s), 'false', 5) | strncmp(s, '0 ', 2)
      res = 0;
    else
      try res = evalin('caller', s); catch
	error(['String value "' s '" cannot be evaluated']);
      end
      try res ~= 0; catch
	error(['String value "' s '" cannot be evaluated reasonably']);
      end
    end
  end
  
% ---------------------------------------------------------------  
% ---------------------------------------------------------------  
function outplot(o)
  persistent lasteval
  acteval = evalin('caller', 'counteval');
  if isempty(lasteval) || acteval < lasteval
    figure(324);
  elseif gcf ~= 324 % prevents repeated raise of figure
    if  ismember(324, findobj('Type', 'figure'))
      set(0, 'CurrentFigure', 324); % geht nur, wenn figure schon exisitiert
    else
      figure(324);
    end
  end
  lasteval = acteval;

  foffset = 1e-99;
  dfit = o.y1(:,1)-min(o.y1(:,1)); 
  dfit(find(dfit<1e-98)) = NaN;
  subplot(2,2,1); hold off; semilogy(o.x,dfit,'-c'); hold on;
  idx = find(o.y1(:,1)>1e-98);  % positive values
  subplot(2,2,1);semilogy(o.x(idx), o.y1(idx,1)+foffset, '.b'); hold on; 
  idx = find(o.y1(:,1) < -1e-98);  % negative values
  subplot(2,2,1);semilogy(o.x(idx), abs(o.y1(idx,1))+foffset,'.r');hold on; 
  subplot(2,2,1);semilogy(o.x,abs(o.y1(:,1))+foffset,'-b'); hold on;
  if size(o.y1, 2) > 3
    %qqq
    subplot(2,2,1);semilogy(o.x,(o.y1(:,4:end)),'-y'); hold on;
    % subplot(2,2,1);semilogy(o.x,(o.y1(:,4:end)),'-m'); hold on;
    % subplot(2,2,1);semilogy(o.x,(o.y1(:,4:end))); hold on;
  end
  subplot(2,2,1);semilogy(o.x,(o.y1(:,3)),'-r'); hold on;
  subplot(2,2,1);semilogy(o.x,(o.y1(:,2)),'-g'); 
  ax = axis;
  text(ax(1), 10^(log10(ax(3))+0.05*(log10(ax(4))-log10(ax(3)))), ...
       [ ' f=' num2str(o.y1(end, 1), '%.15g')]);
  title('abs(f) (blue), f-min(f) (cyan), Sigma (green), Axis Ratio (red)');
  grid on; 

  subplot(2,2,2); hold off; plot(o.x, o.y2,'-'); 
  if size(o.y2, 2) < 100
    ax = axis;
    ax(2) = max(1.05*o.x(end), ax(2));
    axis(ax);
    yy = linspace(ax(3), ax(4), size(o.y2,2))';
    [yyl idx] = sort(o.y2(end,:));
    [muell idx2] = sort(idx);
    hold on;
    plot([o.x(end) ax(2)]', [o.y2(end,:)' yy(idx2)]', '-');
    plot(repmat(o.x(end),2), [ax(3) ax(4)], 'k-');
    for i = 1:length(idx)
      text(ax(2), yy(i), [' ' num2str(idx(i))]);
    end
  end
  title(['Object Variables (' num2str(size(o.y2, 2)) 'D)']);grid on;

  subplot(2,2,3); hold off; semilogy(o.x, o.y3, '-'); 
  if size(o.y2, 2) < 100
    ax = axis; 
    ax(2) = max(1.05*o.x(end), ax(2));
    axis(ax);
    yy = logspace(log10(ax(3)), log10(ax(4)), size(o.y3,2))';
    [yyl idx] = sort(o.y3(end,:));
    [muell idx2] = sort(idx);
    hold on;
    plot([o.x(end) ax(2)]', [o.y3(end,:)' yy(idx2)]', '-');
    plot(repmat(o.x(end),2), [ax(3) ax(4)], 'k-');
    for i = 1:length(idx)
      text(ax(2), yy(i), [' ' num2str(idx(i))]);
    end
  end
  title('Standard Deviations of All Variables'); grid on;
  xlabel('function evaluations'); 

  subplot(2,2,4); semilogy(o.x, o.y4, '-');
  title('Scaling (All Main Axes)'); grid on;
  xlabel('function evaluations'); 
  zoom on; drawnow;
  
% ---------------------------------------------------------------  
% ---------------------------------------------------------------  
% ----- replacements for statistic toolbox functions ------------
% ---------------------------------------------------------------  
% ---------------------------------------------------------------  
function res=myrange(x)
  res = max(x) - min(x);
  
% ---------------------------------------------------------------  
% ---------------------------------------------------------------  
function res = myprctile(inar, perc, idx)
%
% Computes the percentiles in vector perc from vector inar
% returns vector with length(res)==length(perc)
% idx: optional index-array indicating sorted order
%

N = length(inar);
flgtranspose = 0;

% sizes 
if size(perc,1) > 1
  perc = perc';
  flgtranspose = 1;
  if size(perc,1) > 1
    error('perc must not be a matrix');
  end
end
if size(inar, 1) > 1 & size(inar,2) > 1
  error('data inar must not be a matrix');
end
 
% sort inar
if nargin < 3 | isempty(idx)
  [sar idx] = sort(inar);
else
  sar = inar(idx);
end

res = [];
for p = perc
  if p <= 100*(0.5/N)
    res(end+1) = sar(1);
  elseif p >= 100*((N-0.5)/N)
    res(end+1) = sar(N);
  else
    % find largest index smaller than required percentile
    availablepercentiles = 100*((1:N)-0.5)/N;
    i = max(find(p > availablepercentiles));
    % interpolate linearly
    res(end+1) = sar(i) ...
	+ (sar(i+1)-sar(i))*(p - availablepercentiles(i)) ...
	/ (availablepercentiles(i+1) - availablepercentiles(i));

  end
end

if flgtranspose
  res = res';
end


% ---------------------------------------------------------------  
% --------------- OBJECTIVE TEST FUNCTIONS ----------------------  
% ---------------------------------------------------------------  

%%% Unimodal functions

function f=fsphere(x)
  f=sum(x.^2);

function f = fsphereoneax(x)
  f = x(1)^2;
  f = mean(x)^2;
  
function f=frandsphere(x)
  N = size(x,1);
  idx = ceil(N*rand(7,1));
  f=sum(x(idx).^2);

function f=fspherelb0(x, M) % lbound at zero for 1:M needed
  if nargin < 2 M = 0; end
  N = size(x,1);
  % M active bounds, f_i = 1 for x = 0
  f = -M + sum((x(1:M) + 1).^2);
  f = f + sum(x(M+1:N).^2);
  
function f=fspherehull(x)
  % Patton, Dexter, Goodman, Punch
  % in -500..500
  % spherical ridge through zeros(N,1)
  % worst case start point seems x = 2*100*sqrt(N)
  % and small step size
  N = size(x,1);
  f = norm(x) + (norm(x-100*sqrt(N)) - 100*N)^2;
  
function f=fellilb0(x, idxM, scal) % lbound at zero for 1:M needed
  N = size(x,1);
  if nargin < 3 | isempty(scal)
    scal = 100;
  end
  scale=scal.^((0:N-1)/(N-1));
  if nargin < 2 | isempty(idxM)
    idxM = 1:N;
  end
  %scale(N) = 1e0;
  % M active bounds
  xopt = 0.1;
  x(idxM) = x(idxM) + xopt;
  f = scale.^2*x.^2;
  f = f - sum((xopt*scale(idxM)).^2); 
%  f = exp(f) - 1;
%  f = log10(f+1e-19) + 19;

  f = f + 1e-19;
  
function f=fcornersphere(x)
  w = ones(size(x,1));
  w(1) = 2.5; w(2)=2.5;
  idx = x < 0;
  f = sum(x(idx).^2);
  idx = x > 0;
  f = f + 2^2*sum(w(idx).*x(idx).^2);
  
function f=fsectorsphere(x, scal)
%
% This is deceptive for cumulative sigma control in large dimension:
% The strategy (initially) diverges for N=50 and popsize = 150.  (Even
% for cs==1 this can be observed for larger settings of N and
% popsize.) The reason is obvious from the function topology. 
% Divergence can be avoided by setting boundaries or adding a
% penalty for large ||x||. Then, convergence can be observed again. 
% Conclusion: for popsize>N cumulative sigma control is not completely
% reasonable, but I do not know better alternatives.
%
  if nargin < 2 | isempty (scal)
    scal = 1e3;
  end
  f=sum(x.^2);
  idx = find(x<0);
  f = f + (scal-1)^2 * sum(x(idx).^2);
  
function f=fstepsphere(x, scal)
  if nargin < 2 | isempty (scal)
    scal = 1e0;
  end
  N = size(x,1);
  f=1e-11+sum(scal.^((0:N-1)/(N-1))*floor(x+0.5).^2);
  f=1e-11+sum(floor(scal.^((0:N-1)/(N-1))'.*x+0.5).^2);
%  f=1e-11+sum(floor(x+0.5).^2);

function f=fstep(x)
  % in -5.12..5.12 (bounded)
  N = size(x,1);
  f=1e-11+6*N+sum(floor(x));

function f=flnorm(x, scal, e)
if nargin < 2 | isempty(scal)
  scal = 1;
end
if nargin < 3 | isempty(e)
  e = 1;
end
if e==inf
  f = max(abs(x));
else
  N = size(x,1);
  scale = scal.^((0:N-1)/(N-1))';
  f=sum(abs(scale.*x).^e);
end

function f=fneumaier3(x) 
  % in -n^2..n^2
  % x^*-i = i(n+1-i)
  N = size(x,1);
%  f = N*(N+4)*(N-1)/6 + sum((x-1).^2) - sum(x(1:N-1).*x(2:N));
  f = sum((x-1).^2) - sum(x(1:N-1).*x(2:N));
  
function f=fchangingsphere(x)
  N = size(x,1);
  global scale_G; global count_G; if isempty(count_G) count_G=-1; end
  count_G = count_G+1;
  if mod(count_G,10) == 0
    scale_G = 10.^(2*rand(1,N));
  end
  %disp(scale(1));
  f = scale_G*x.^2;
  
function f= flogsphere(x)
 f = 1-exp(-sum(x.^2));
  
function f= fexpsphere(x)
 f = exp(sum(x.^2)) - 1;
  
function f=fbaluja(x)
  % in [-0.16 0.16]
  y = x(1);
  for i = 2:length(x)
    y(i) = x(i) + y(i-1);
  end
  f = 1e5 - 1/(1e-5 + sum(abs(y)));

function f=fschwefel(x)
  f = 0;
  for i = 1:size(x,1),
    f = f+sum(x(1:i))^2;
  end

function f=fcigar(x)
  f = x(1)^2 + 1e6*sum(x(2:end).^2);
  
function f=fcigtab(x)
  f = x(1)^2 + 1e8*x(end)^2 + 1e4*sum(x(2:(end-1)).^2);
  
function f=ftablet(x)
  f = 1e6*x(1)^2 + sum(x(2:end).^2);
  
function f=felli(x, lgscal, expon, expon2)
  % lgscal: log10(axis ratio)
  % expon: x_i^expon, sphere==2
  N = size(x,1); if N < 2 error('dimension must be greater one'); end
  if nargin < 2 | isempty(lgscal), lgscal = 3; end
  if nargin < 3 | isempty(expon), expon = 2; end
  if nargin < 4 | isempty(expon2), expon2 = 1; end

  f=((10^(lgscal*expon)).^((0:N-1)/(N-1)) * abs(x).^expon)^(1/expon2);
%  if rand(1,1) > 0.015
%    f = NaN;
%  end

function f=fellii(x, scal)
  N = size(x,1); if N < 2 error('dimension must be greater one'); end
  if nargin < 2
    scal = 1;
  end
  f= (scal*(1:N)).^2 * (x).^2;

function f=fellirot(x)
  N = size(x,1);
  global ORTHOGONALCOORSYSTEM_G
  if isempty(ORTHOGONALCOORSYSTEM_G) ...
	| isempty(ORTHOGONALCOORSYSTEM_G{N})
    coordinatesystem(N);
  end
  f = felli(ORTHOGONALCOORSYSTEM_G{N}*x);
  
function coordinatesystem(N)
  if nargin < 1 | isempty(N)
    arN = 2:30;
  else
    arN = N;
  end
  global ORTHOGONALCOORSYSTEM_G
  ORTHOGONALCOORSYSTEM_G{1} = 1; 
  for N = arN
    ar = randn(N,N);
    for i = 1:N 
      for j = 1:i-1
	ar(:,i) = ar(:,i) - ar(:,i)'*ar(:,j) * ar(:,j);
      end
      ar(:,i) = ar(:,i) / norm(ar(:,i));
    end
    ORTHOGONALCOORSYSTEM_G{N} = ar; 
  end

function f=fplane(x)
  f=x(1);

function f=ftwoaxes(x)
  f = sum(x(1:floor(end/2)).^2) + 1e6*sum(x(floor(1+end/2):end).^2);

function f=fparabR(x)
  f = -x(1) + 100*sum(x(2:end).^2);

function f=fsharpR(x)
  f = abs(-x(1)) + 30*norm(x(2:end));
  
function f=frosen(x)
  if size(x,1) < 2 error('dimension must be greater one'); end
  f = 100*sum((x(1:end-1).^2 - x(2:end)).^2) + sum((x(1:end-1)-1).^2);
  
function f=frosenmodif(x)
  f = 74 + 100*(x(2)-x(1)^2)^2 + (1-x(1))^2 ...
      - 400*exp(-sum((x+1).^2)/2/0.05);
  
function f=fschwefelrosen1(x)
  % in [-10 10] 
  f=sum((x.^2-x(1)).^2 + (x-1).^2);
  
function f=fschwefelrosen2(x)
  % in [-10 10] 
  f=sum((x(2:end).^2-x(1)).^2 + (x(2:end)-1).^2);

function f=fdiffpow(x)
  N = size(x,1); if N < 2 error('dimension must be greater one'); end
  f=sum(abs(x).^(2+10*(0:N-1)'/(N-1)));

%%% Multimodal functions 

function f=fackley(x)
  % -32.768..32.768
  % Adding a penalty outside the interval is recommended,  
  % because for large step sizes, fackley imposes like frand
  % 
  N = size(x,1); 
  f = 20-20*exp(-0.2*sqrt(sum(x.^2)/N)); 
  f = f + (exp(1) - exp(sum(cos(2*pi*x))/N));
  % add penalty outside the search interval
  f = f + sum((x(x>32.768)-32.768).^2) + sum((x(x<-32.768)+32.768).^2);
  
function f = fbohachevsky(x)
 % -15..15
  f = sum(x(1:end-1).^2 + 2 * x(2:end).^2 - 0.3 * cos(3*pi*x(1:end-1)) ...
	  - 0.4 * cos(4*pi*x(2:end)) + 0.7);
  
function f=fconcentric(x)
  % in  +-600
  s = sum(x.^2);
  f = s^0.25 * (sin(50*s^0.1)^2 + 1);

function f=fgriewank(x)
  % in [-600 600]
  N = size(x,1);
  f = 1 - prod(cos(x'./sqrt(1:N))) + sum(x.^2)/4e3;
  % f = f + 1e4*sum(x(abs(x)>5).^2);
  % if sum(x(abs(x)>5).^2) > 0
  %   f = 1e4 * sum(x(abs(x)>5).^2) + 1e8 * sum(x(x>5)).^2;
  % end
  
function f=frastrigin(x, scal, skewfac, skewstart, amplitude)
  if nargin < 5 | isempty(amplitude)
    amplitude = 10;
  end
  if nargin < 4 | isempty(skewstart)
    skewstart = 0;
  end
  if nargin < 3 | isempty(skewfac)
    skewfac = 1;
  end
  if nargin < 2 | isempty(scal)
    scal = 1;
  end
  N = size(x,1); 
  scale = 1;
  if N > 1
    scale=scal.^((0:N-1)'/(N-1));
  end
  % simple version: 
  % f = amplitude*(N - sum(cos(2*pi*(scale.*x)))) + sum((scale.*x).^2);

  % skew version: 
  y = scale.*x;
  idx = find(x > skewstart);
  if ~isempty(idx)
    y(idx) =  skewfac*x(idx);
  end
  f = amplitude * (N-sum(cos(2*pi*y))) + sum(y.^2);
  
function f = fschaffer(x)
 % -100..100
  N = size(x,1);
  s = x(1:N-1).^2 + x(2:N).^2;
  f = sum(s.^0.25 .* (sin(50*s.^0.1).^2+1));

function f=fschwefelmult(x)
  % -500..500
  N = size(x,1); 
  f = 418.9829*N - 1.27275661e-5*N - sum(x.*sin(sqrt(abs(x))));
  f = f + sum(x(abs(x)>500).^2);
  
function f=ftwomax(x)
  % Boundaries at +/-5
  N = size(x,1); 
  f = -abs(sum(x)) + 5*N;

function f=ftwomaxtwo(x)
  % Boundaries at +/-10
  N = size(x,1); 
  f = abs(sum(x));
  if f > 30
    f = f - 30;
  end
  f = -f;

function f=frand(x)
  f=1/(1-rand) - 1;

% Changes: 
% 05/09: Raise of figure and waiting for first plots improved
% 05/01: Function coordinatesystem cleaned up. 
% 05/01: Function prctile, which requires the statistics toolbox,
%        replaced by myprctile. 
% 05/01: Option warnonequalfunctionvalues included. 
% 04/12: Decrease of sigma removed. Problems on fsectorsphere can 
%        be addressed better by adding search space boundaries. 
% 04/12: Boundary handling simpyfied. 
% 04/12: Bug when stopping criteria tolx or tolupx are vectors. 
% 04/11: Three input parameters are obligatory now. 
% 04/11: Bug in boundary handling removed: Boundary weights can decrease now. 
% 04/11: Normalization for boundary weights scale changed. 
% 04/11: VerboseModulo option bug removed. Documentation improved. 
% 04/11: Condition for increasing boundary weights changed.
% 04/10: Decrease of sigma when fitness is getting consistenly
%        worse. Addresses the problems appearing on fsectorsphere for
%        large population size.
% 04/10: VerboseModulo option included. 
% 04/10: Bug for condition for increasing boundary weights removed.
% 04/07: tolx depends on initial sigma to achieve scale invariance
%        for this stopping criterion. 
% 04/06: Objective function value NaN is not counted as function
%        evaluation and invokes resampling of the search point. 
% 04/06: Error handling for eigenvalue beeing zero (never happens
%        with default parameter setting)
% 04/05: damps further tuned for large mueff 
%      o Details for stall of pc-adaptation added (variable hsig 
%        introduced). 
% 04/05: Bug in boundary handling removed: A large initial SIGMA was
%        corrected not until *after* the first iteration, which could
%        lead to a complete failure.
% 04/05: Call of function range (works with stats toolbox only) 
%        changed to myrange. 
% 04/04: Parameter cs depends on mueff now and damps \propto sqrt(mueff)
%        instead of \propto mueff. 
%      o Initial stall to adapt C (flginiphase) is removed and
%        adaptation of pc is stalled for large norm(ps) instead.
%      o Returned default options include documentation. 
%      o Resume part reorganized.
% 04/03: Stopflag becomes cell-array. 

% ---------------------------------------------------------------
% CMA-ES: Evolution Strategy with Covariance Matrix Adaptation for
% nonlinear function minimization. To be used under the terms of the
% GNU General Public License (http://www.gnu.org/copyleft/gpl.html).
% Author: Nikolaus Hansen, 2001/3. e-mail: hansen@bionik.tu-berlin.de
% URL:http://www.bionik.tu-berlin.de/user/niko
% References: See below. 
% ---------------------------------------------------------------
%
% GENERAL PURPOSE: The CMA-ES (Evolution Strategy with Covariance
% Matrix Adaptation) is a robust search method which should be
% applied, if derivative based methods, e.g. quasi-Newton BFGS or
% conjucate gradient, (supposably) fail due to a rugged search
% landscape (e.g. noise, local optima, outlier, etc.). On smooth
% landscapes CMA-ES is roughly ten times slower than BFGS. For up to
% N=10 variables even the simplex direct search method (Nelder & Mead)
% is often faster, but far less robust than CMA-ES.  To see the
% advantage of the CMA, it will usually take at least 30*N and up to
% 300*N function evaluations, where N is the search problem dimension.
% On considerably hard problems the complete search (a single run) is
% expected to take at least 30*N^2 and up to 300*N^2 function
% evaluations.
%
% SOME MORE COMMENTS: 
% The adaptation of the covariance matrix (e.g. by the CMA) is
% equivalent to a general linear transformation of the problem
% coding. Nevertheless every problem specific knowlegde about the best
% linear transformation should be exploited before starting the
% search. That is, an appropriate a priori transformation should be
% applied to the problem. This also makes the identity matrix as
% initial covariance matrix the best choice.
%
% The strategy parameter lambda (population size, opts.PopSize) is the
% preferred strategy parameter to play with.  If results with the
% default strategy are not satisfactory, increase the population
% size. (Remark that the crucial parameter mu (opts.ParentNumber) is
% increased proportionally to lambda). This will improve the
% strategies capability of handling noise and local minima. We
% recomment successively increasing lambda by a factor of about three,
% starting with initial values between 5 and 20. Casually, population
% sizes even beyond 1000+100*N can be sensible.
%
%
% ---------------------------------------------------------------
%%% REFERENCES
%
% The equation numbers refer to 
% Hansen, N. and S. Kern (2004). Evaluating the CMA Evolution
% Strategy on Multimodal Test Functions.  Eighth International
% Conference on Parallel Problem Solving from Nature PPSN VIII,
% Proceedings, pp. 282-291, Berlin: Springer. 
% (http://www.bionik.tu-berlin.de/user/niko/ppsn2004hansenkern.pdf)
% 
% Further references:
% Hansen, N. and A. Ostermeier (2001). Completely Derandomized
% Self-Adaptation in Evolution Strategies. Evolutionary Computation,
% 9(2), pp. 159-195.
% (http://www.bionik.tu-berlin.de/user/niko/cmaartic.pdf).
%
% Hansen, N., S.D. Mueller and P. Koumoutsakos (2003). Reducing the
% Time Complexity of the Derandomized Evolution Strategy with
% Covariance Matrix Adaptation (CMA-ES). Evolutionary Computation,
% 11(1).  (http://mitpress.mit.edu/journals/pdf/evco_11_1_1_0.pdf).
%

