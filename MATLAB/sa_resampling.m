%  The code is a Matlab version by Martin M. Andreasen of the code 
%  written by Alan Miller and used in "Global Optimization of Statistical
%  Functions with Simulated Annealing," Goffe, Ferrier and
%  Rogers, Journal of Econometrics, vol. 60, no. 1/2, Jan./Feb. 1994, pp.
%  65-100. To make Simulated Annealing work for objective functions defined on a
%  non-convex domain, Martin M. Andreasen have extended Simulated Annealing
%  with resampling and jumps  
%  Latest version: The spring of 2007

function [x,t,vm,xopt,fopt,nacc,nfcnev,nobds,ier] = sa_resampling(fitfun,n,x,max,rt,eps,ns,nt,neps,maxevl,...
    maxresample,maxNaNJump,lb,ub,c,iprint,t,vm,Save_To_File);

%*************************************************************************
%  EXAMPLE: How to use this minimizer
%  clc;                       % Clear the window
%  n = 4;                     % The objective function has 4 parameters
%  x0 = [80 -10 10 20];       % The starting values
%  max = 0;                   % We do minimization
%  rt = 0.85;                 % The reduction rate in temperature
%  eps = 1e-6;                % eps
%  ns = 10;                   % Number of cycles
%  nt = 20;                   % Number of random walkers
%  neps = 4;                  % Stopping criteria
%  maxevl = 1e5;              % Maximum number of function evaluations
%  maxresample = 15;          % Maximum number of times we resample a new step size
%  maxNaNJump = 5;            % Maximum number of times we jump for a given step size
%  lb = -[100 100 100 100];   % Lower bound
%  ub =  [100 100 100 100];   % Upper bound
%  c = [2 2 2 2];             % Default setting for adjusting vm - the stepsizes
%  iprint = 1;                % Printing property
%  t = 100;                   % Initial temperature
%  vm = [10 5 5 5] ;          % Initial stepsizes

% [x,t,vm,xopt,fopt,nacc,nfcnev,nobds,ier] = sa_resampling(@test_function,n,x0,max,rt,eps,ns,nt,neps,maxevl,   maxresample,maxNaNJump,lb,ub,c,iprint,t,vm);

%  THE OBJECTIVE FUNCTION
%  function f=test_function(x)
%  f = 100*(0-x(1))^2 + 20*(0-x(2))^2+(0-x(3))^2+(0-x(4))^4;
%  if  65 < x(1) & x(1) < 75
%      f = NaN;
%  end

%*************************************************************************

%  Synopsis:
%  A very quick (perhaps too quick) overview of SA:
%     SA tries to find the global optimum of an N dimensional function.
%  It moves both up and downhill and as the optimization process
%  proceeds, it focuses on the most promising area.
%     To start, it randomly chooses a trial point within the step length
%  VM (a vector of length N) of the user selected starting point. The
%  function is evaluated at this trial point and its value is compared
%  to its value at the initial point.
%     In a maximization problem, all uphill moves are accepted and the
%  algorithm continues from that trial point. Downhill moves may be
%  accepted; the decision is made by the Metropolis criteria. It uses T
%  (temperature) and the size of the downhill move in a probabilistic
%  manner. The smaller T and the size of the downhill move are, the more
%  likely that move will be accepted. If the trial is accepted, the
%  algorithm moves on from that point. If it is rejected, another point
%  is chosen instead for a trial evaluation.
%     Each element of VM periodically adjusted so that half of all
%  function evaluations in that direction are accepted.
%     A fall in T is imposed upon the system with the RT variable by
%  T(i+1) = RT*T(i) where i is the ith iteration. Thus, as T declines,
%  downhill moves are less likely to be accepted and the percentage of
%  rejections rise. Given the scheme for the selection for VM, VM falls.
%  Thus, as T declines, VM falls and SA focuses upon the most promising
%  area for optimization.

%  The importance of the parameter T:
%     The parameter T is crucial in using SA successfully. It influences
%  VM, the step length over which the algorithm searches for optima. For
%  a small intial T, the step length may be too small; thus not enough
%  of the function might be evaluated to find the global optima. The user
%  should carefully examine VM in the intermediate output (set IPRINT =
%  1) to make sure that VM is appropriate. The relationship between the
%  initial temperature and the resulting step length is function
%  dependent.
%     To determine the starting temperature that is consistent with
%  optimizing a function, it is worthwhile to run a trial run first. Set
%  RT = 1.5 and T = 1.0. With RT > 1.0, the temperature increases and VM
%  rises as well. Then select the T that produces a large enough VM.

%  This version of Simulated Annealing is developed by Corana et al.'s article "Minimizing
%  Multimodal Functions of Continuous Variables with the "Simulated Annealing"
%  Algorithm" in the September 1987 (vol. 13, no. 3, pp. 262-280) issue of
%  the ACM Transactions on Mathematical Software. For modifications to the 
%  algorithm and many details on its use, (particularly for econometric applications) 
%  see Goffe, Ferrier and Rogers, "Global Optimization of Statistical Functions with
%  Simulated Annealing," Journal of Econometrics, vol. 60, no. 1/2,
%  Jan./Feb. 1994, pp. 65-100.
%  For more information, contact
%              Bill Goffe
%              Department of Economics and International Business
%              University of Southern Mississippi
%              Hattiesburg, MS  39506-5072
%              (601) 266-4484 (office)
%              (601) 266-4920 (fax)
%              bgoffe@whale.st.usm.edu (Internet)
%  Martin M. Andreasen have further extended the algorithm with resampling and jumps
%  in case the objective function can not be evaluated in a given point. 
%  Notice, if undefined points for the objective function are given a very 
%  large value then the resampling and jump procedure  never starts.
%  I.e. this code nests the version of Simulated Annealing by Goffe, Ferrier
%  and Rogers

%  As far as possible, the parameters here have the same name as in
%  the description of the algorithm on pp. 266-8 of Corana et al.
%  In this description, SP is single precision, DP is double precision,
%  INT is integer, L is logical and (N) denotes an array of length n.
%  Thus, DP(N) denotes a double precision array of length n.

%  Input Parameters:
%    Note: The suggested values generally come from Corana et al. To
%          drastically reduce runtime, see Goffe et al., pp. 90-1 for
%          suggestions on choosing the appropriate RT and NT.
%    N - Number of variables in the function to be optimized. (INT)
%    X - The starting values for the variables of the function to be
%        optimized. (DP(N))
%    MAX - Denotes whether the function should be maximized or minimized.
%          A true value denotes maximization while a false value denotes
%          minimization.  Intermediate output (see IPRINT) takes this into
%          account. (L)
%    RT - The temperature reduction factor.  The value suggested by
%         Corana et al. is .85. See Goffe et al. for more advice. (DP)
%    EPS - Error tolerance for termination. If the final function
%          values from the last neps temperatures differ from the
%          corresponding value at the current temperature by less than
%          EPS and the final function value at the current temperature
%          differs from the current optimal function value by less than
%          EPS, execution terminates and IER = 0 is returned. (EP)
%    NS - Number of cycles.  After NS*N function evaluations, each element of
%         VM is adjusted so that approximately half of all function evaluations
%         are accepted.  The suggested value is 20. (INT)
%    NT - Number of iterations before temperature reduction. After
%         NT*NS*N function evaluations, temperature (T) is changed
%         by the factor RT.  Value suggested by Corana et al. is
%         MAX(100, 5*N).  See Goffe et al. for further advice. (INT)
%    NEPS - Number of final function values used to decide upon termi-
%           nation.  See EPS.  Suggested value is 4. (INT)
%    MAXEVL - The maximum number of function evaluations.  If it is
%             exceeded, IER = 1. (INT)
%    MAXRESAMPLE - The maximum number of times we resample in each move to get
%                a point where the object function is defined.(INT) 
%    MAXNANJUMP  - The maximum number of times we scale the innovation in x(i)
%                  (i.e. jump) to find a defined point
%    LB - The lower bound for the allowable solution variables. (DP(N))
%    UB - The upper bound for the allowable solution variables. (DP(N))
%         If the algorithm chooses X(I) .LT. LB(I) or X(I) .GT. UB(I),
%         I = 1, N, a point is from inside is randomly selected. This
%         This focuses the algorithm on the region inside UB and LB.
%         Unless the user wishes to concentrate the search to a particular
%         region, UB and LB should be set to very large positive
%         and negative values, respectively.  Note that the starting
%         vector X should be inside this region.  Also note that LB and
%         UB are fixed in position, while VM is centered on the last
%         accepted trial set of variables that optimizes the function.
%    C - Vector that controls the step length adjustment.  The suggested
%        value for all elements is 2.0. (DP(N))
%    IPRINT - controls printing inside SA. (INT)
%             Values: 0 - Only termination information is printed.
%                     1 - Function value for the starting value and
%                         summary results before each temperature
%                         reduction. This includes the optimal
%                         function value found so far, the total
%                         number of moves (broken up into uphill,
%                         downhill, accepted and rejected), the
%                         number of out of bounds trials, the
%                         number of new optima found at this
%                         temperature, the current optimal X and
%                         the step length VM. Note that there are
%                         N*NS*NT function evalutations before each
%                         temperature reduction. Finally, notice is
%                         is also given upon achieveing the termination
%                         criteria.
%                     2 - Each new step length (VM), the current optimal
%                         X (XOPT) and the current trial X (X). This
%                         gives the user some idea about how far X
%                         strays from XOPT as well as how VM is adapting
%                         to the function.
%                     3 - Each function evaluation, its acceptance or
%                         rejection and new optima. For many problems,
%                         this option will likely require a small tree
%                         if hard copy is used. This option is best
%                         used to learn about the algorithm. A small
%                         value for MAXEVL is thus recommended when
%                         using IPRINT = 4.
%             Suggested value: 1
%             Note: For a given value of IPRINT, the lower valued
%                   options (other than 0) are utilized.
%    ISEED1 - The first seed for the random number generator RANMAR.
%             0 <= ISEED1 <= 31328. (INT)
%    ISEED2 - The second seed for the random number generator RANMAR.
%             0 <= ISEED2 <= 30081. Different values for ISEED1
%             and ISEED2 will lead to an entirely different sequence
%             of trial points and decisions on downhill moves (when
%             maximizing).  See Goffe et al. on how this can be used
%             to test the results of SA. (INT)

%  Input/Output Parameters:
%    T - On input, the initial temperature. See Goffe et al. for advice.
%        On output, the final temperature. (DP)
%    VM - The step length vector. On input it should encompass the region of
%         interest given the starting value X.  For point X(I), the next
%         trial point is selected is from X(I) - VM(I)  to  X(I) + VM(I).
%         Since VM is adjusted so that about half of all points are accepted,
%         the input value is not very important (i.e. is the value is off,
%         SA adjusts VM to the correct value). (DP(N))

%  Output Parameters:
%    XOPT - The variables that optimize the function. (DP(N))
%    FOPT - The optimal value of the function. (DP)
%    NACC - The number of accepted function evaluations. (INT)
%    NFCNEV - The total number of function evaluations. In a minor
%             point, note that the first evaluation is not used in the
%             core of the algorithm; it simply initializes the
%             algorithm. (INT).
%    NOBDS - The total number of trial function evaluations that
%            would have been out of bounds of LB and UB. Note that
%            a trial point is randomly selected between LB and UB. (INT)
%    IER - The error return number. (INT)
%          Values: 0 - Normal return; termination criteria achieved.
%                  1 - Number of function evaluations (NFCNEV) is
%                      greater than the maximum number (MAXEVL).
%                  2 - The starting value (X) is not inside the
%                      bounds (LB and UB).
%                  3 - The initial temperature is not positive.
%                  99 - Should not be seen; only used internally.

%  Required Subroutines (included):
%    PRTVEC - Prints vectors.
%    PRT1 ... PRT10 - Prints intermediate output.
%    FCN - Function to be optimized. The form is
%            function FCN(N, X, F, Error_mes)




%Clearing the window for reporting
clc;

%Initialize the random number generator.
rand('state', sum(100*clock));

%Set initial values
nacc      = 0;           %Total number of accepted function evaluations.
nobds     = 0;           %Total number of function evaluations that would have been outside (lb,ub)
nfcnev    = 0;           %Total number of function evaluations
ier       = 99;
xopt(1:n) = x(1:n);
nacp      = zeros(1,n);
fstar     = 1.0E+20*ones(1,neps);

%If the initial temperature is not positive, notify the user and return to the calling routine.
if t <= 0.0 
    display(' THE INITIAL TEMPERATURE IS NOT POSITIVE. RESET THE VALUE OF t  ');
    ier = 3;
    xopt=NaN;
    fopt=NaN;
    return;
end

%If the initial value is out of bounds, notify the user and return to the calling routine.
for i=1:n
    if x(i) > ub(i) | x(i) < lb(i) 
        disp(['x(i)= ' num2str(x(i)) ' PARAMETER NUMBER ' num2str(i) ' ub(i) ' num2str(ub(i)) ' lb(i) ' num2str(lb(i)) ]);
        %CALL prt1()
        disp('THE STARTING VALUE (X) IS OUTSIDE THE BOUNDS (lb AND ub)');
        disp('Execution terminated without any optimization.');
        disp('Respecify x, ub OR lb so that lb(i) < x(i) < ub(i), i = 1, n');
        ier = 2;
        xopt=NaN;
        fopt=NaN;
        return;
    end
end

%Evaluate the function with input x and return value as f.
f = feval(fitfun,x);
% if f = NaN then give f the value 1e20 
if isnan(f) == 1
    f = 1.E20;
end

%If the function is to be minimized, switch the sign of the function.
%Note that all intermediate and final output switches the sign back
%to eliminate any possible confusion for the user.
if max == 1
    f = f;
else
    f = -f;
end
nfcnev = nfcnev + 1;
fopt = f;
fstar(1) = f;

%The inital conditions for the routine are displayed
if iprint >= 1
   disp('****  SIMULATED ANNEALING  ****');
   disp(' ');
   disp(['NUMBER OF PARAMETERS                   ', num2str(n)]);
   if max == 1
       disp('OBJECTIVE FUNCTION IS BEING            MAXIMIZED');
   else
       disp('OBJECTIVE FUNCTION IS BEING            MINIMIZED');
   end
   disp(['INITIAL TEMPERATURE                    ', num2str(t)]);
   disp(['REDUCTION RATE IN TEMPERATURE          ', num2str(rt)]);
   disp(['EPS                                    ', num2str(eps)]);
   disp(['NUMBER OF neps                         ', num2str(neps)]);
   disp(['NUMBER OF CYCLES (ns)                  ', num2str(ns)]);
   disp(['NUMBER OF RANDOM WALKERS (nt)          ', num2str(nt)]);
   disp(['MAXRESAMPLE                            ', num2str(maxresample)]);
   disp(['MAXNANJUMP                             ', num2str(maxNaNJump)]);
   disp(['MAXIMUM NUMBER OF FUNCTION EVALUATIONS ', num2str(maxevl)]);
   disp(['DISPLAY PROPERTIES                     ', num2str(iprint)]);
   disp(' ');
   prtvec(x, n, 'STARTING VALUES');
   prtvec(vm, n, 'INITIAL STEP LENGTH');
   prtvec(lb, n, 'LOWER BOUND');
   prtvec(ub, n, 'UPPER BOUND');
   prtvec(c, n, 'C VECTOR');
   disp('****   BEFORE STARTING SIMULATED ANNEALING   ****');
end


%Start the main loop. Note that it terminates if (i) the algorithm
%succesfully optimizes the function or (ii) there are too many
%function evaluations (more than MAXEVL).
quit    = 0;
mark_10 = 0;

while quit == 0
nup         = 0;          %Number of accepted function evaluations going 'up' for a given temperature
index       = 0;          %Index for the number of times we draw a new xp for each i=1,...,n
nrej        = 0;          %Number of rejected function evalutions for a given temperature
nnew        = 0;          %Number of new optimums for a given temperature
ndown       = 0;          %Number of accepted function evalutions going 'down' for a given temperature
lnobds      = 0;          %Number of function evaluations that would have been outside (lb,ub) for a given temperature
num_no_move = zeros(1,n); %Number of times where we cannot find a move, in direction i=1,...,n


for m=1:nt
   for j=1:ns
       for h=1:n
          mark_10 = 0;
          mark_40 = 0;
          %Generate XP, the trial value of X. Note use of VM to choose XP.
          while mark_10 == 0   
             for i=1:n
                if i == h
                   x_step = (rand*2 - 1)*vm(i);
                   xp(i) = x(i) + x_step;
                else
                   xp(i) = x(i);
                end
        
                %If XP is out of bounds, select a point in bounds for the trial.
                if xp(i) < lb(i) | xp(i) > ub(i) 
                   xp(i) = lb(i) + (ub(i) - lb(i))*rand;
                   lnobds = lnobds + 1;
                   nobds = nobds + 1;
                   if iprint >= 3
                      prt3(max, n, xp, x, f); 
                   end
                end 
             end % loop index by i
      
             %Evaluate the function with the trial point XP and return as FP.
             fp = feval(fitfun, xp);
             if max == 1
                fp = fp;
             else
                fp = -fp;
             end
             nfcnev = nfcnev + 1;
      
             %Checking for errors when evaluating the function at xp. If no error we proceede as
             %in the normal algoritme. 
             if isnan(fp) == 0
                if iprint >= 3
                   disp('NO NEED FOR RESAMPLING');
                   prt4(max, n, xp, x, fp, f);
                end
                mark_10 = 1; %We leave the while loop controlled by mark_10
                
                % If too many function evaluations occur, terminate the algorithm.
                if nfcnev >= maxevl 
                   %CALL prt5()
                   disp(' ');
                   disp('TOO MANY FUNCTION EVALUATIONS; CONSIDER INCREASING maxevl OR eps, OR');
                   disp('DECREASING nt OR rt. THESE RESULTS ARE LIKELY TO BE POOR');
                   disp(' ');
                   if max == 1
                      f = f;
                   else
                      f = -f;
                   end
                   ier = 1;
                   quit = 1;
                   if iprint >= 1
                      disp('****   RESULTS AFTER SIMULATED ANNEALING   ****');
                      prtvec(xopt, n, 'THE BEST POINT');
                      prtvec(vm, n, 'FINAL STEP LENGTH');
                      disp(['OPTIMAL FUNCTION VALUE:             ', num2str(fopt)]);
                      disp(['NUMBER OF FUNCTION EVALUATIONS:     ', num2str(nfcnev)]);
                      disp(['NUMBER OF ACCEPTED EVALUATIONS:     ', num2str(nacc)]);
                      disp(['NUMBER OF OUT OF BOUND EVALUATIONS: ', num2str(nobds)]);
                      disp(['FINAL TEMP:                         ', num2str(t)]);
                      disp(['IER:                                ', num2str(ier)]);
                   end 
                   return;
                end

             %We start the resampling procedure with jumps to get a xp without error, 
             else
                if iprint >=3 & index == 0
                   disp('RESAMPLING STARTS'); 
                end
                index = index + 1;
                if iprint >= 3 & index <= maxresample
                   disp(['NUMBER OF TIMES WE HAVE RESAMPLED THE STEP SIZE ', num2str(index)]); 
                end
                if index <= maxresample
                   %We look for a point without error by scaling the step in the h'te direction
                   indexNaN = 0;
                   while indexNaN < maxNaNJump
                      indexNaN = indexNaN + 1;
                      xp(h) = x(h)+x_step*(1+indexNaN);
                      if iprint >= 3
                          disp(['JUMP NUMBER ', num2str(indexNaN)]);
                      end
                      %Controlling that we do not jump outside lb and ub
                      if xp(h) < lb(h) | xp(h) > ub(h) 
                         lnobds = lnobds + 1;
                         nobds = nobds + 1;
                         if iprint >= 3 
                            prt3(max, n, xp, x, f);
                            disp('JUMP OUTSIDE LB and UB. DRAW A NEW STEP VALUE');
                         end
                         indexNaN = maxNaNJump + 1; %We leave the while loop indexed by indexNaN
                      end

                      if indexNaN < maxNaNJump +1;  %Only enter this if-statement if lb(h)<= xp(h) <= (ub)
                         fp = feval(fitfun, xp);
                         if max == 1
                            fp = fp;
                         else
                            fp = -fp;
                         end
                         nfcnev = nfcnev + 1;
                         if iprint >= 3 
                            prt4(max, n, xp, x, fp, f);
                         end
               
                         if isnan(fp) == 0
                            %We have found a point without error
                            if iprint >= 3
                                disp('FOUND A POINT WITHOUT ERROR. END OF RESAMPLING');
                            end
                            indexNaN = maxNaNJump + 1; %We leave the while loop indexed by indexNaN
                            mark_10 = 1;               %We leave the while loop controlled by mark_10
                         end
                      end
                   end %For while loop indexed by indexNaN
                else 
                   %We cannot find a point without error. Hence we move to the next parameter or (if h=n) start on a new cycle
                   num_no_move(h) = num_no_move(h) + 1;
                   mark_10 = 1;     %We leave the while loop indexed by mark_10
                   mark_40 = 1;     %We cannot find a point without error
                   if iprint >=3 
                       disp(' ');
                       disp('NO POINTS FOUND WITHOUT ERROR. END OF RESAMPLING.');
                       disp('GO TO NEXT DIRECTION OR A NEW CYCLE');
                       disp(' ');
                   end
                end %for the if-statment: if index <= maxresample
             end %for the big if-statement: if isnan(fp) == 0
          end %for while statement indexed by mark_10
          if mark_40 == 1
             %Increase h by 1 if possible, else go to the next cycle
             index   = 0;
             continue; 
          end

          %Resetting the index for the number of times we draw a new xp
          index = 0;
      
          %Accept the new point if the function value increases - for maximum. Notice, if 
          %minimizing then we are maximizing -f
          if fp >= f 
             if iprint >= 3
                disp('POINT ACCEPTED');
                disp(' ');
                disp(['NUMBER OF FUNCTION EVALUATIONS ', num2str(nfcnev)]);
                disp(['RANDOM WARKER NUMBER ', num2str(m), ' IN CYCLE NUMBER ', num2str(j) ,' FOR PARAMETER NUMBER ', num2str(h)]);
             end
             x(1:n) = xp(1:n);
             f = fp;
             nacc = nacc + 1;
             nacp(h) = nacp(h) + 1;
             nup = nup + 1;
        
             %If greater than any other point, record as new optimum.
             if fp > fopt
                if iprint >= 3
                   disp('NEW OPTIMUM');
                end
                xopt(1:n) = xp(1:n);
                fopt = fp;
                nnew = nnew + 1;
             end
             %if the point is lower, use the Metropolis criteria to decide on acceptance or rejection.
          else
             p = exp((fp - f)/t);
             pp = rand;
             if pp < p
                if iprint >= 3 
                   prt6(max);
                   disp(' ');
                   disp(['NUMBER OF FUNCTION EVALUATIONS ', num2str(nfcnev)]);
                   disp(['RANDOM WARKER NUMBER ', num2str(m), ' IN CYCLE NUMBER ',num2str(j),' FOR PARAMETER NUMBER ',num2str(h)]); 
                end
                x(1:n) = xp(1:n);
                f = fp;
                nacc = nacc + 1;
                nacp(h) = nacp(h) + 1;
                ndown = ndown + 1;
             else
                nrej = nrej + 1;
                if iprint >= 3 
                   prt7(max);
                   disp(' ');
                   disp(['NUMBER OF FUNCTION EVALUATIONS ', num2str(nfcnev)]);
                   disp(['RANDOM WARKER NUMBER ', num2str(m), ' IN CYCLE NUMBER ',num2str(j),' FOR PARAMETER NUMBER ',num2str(h)]); 
                end
             end
          end %for if-statement related to the Metropolis Criteria
       end %for the loop indexed by h 
   end %for the loop indexed by j - the number of cycles
   
   %Adjust VM so that approximately half of all evaluations are accepted.
   for i = 1:n
      ratio = nacp(i)/ns;
      if ratio > 0.6
         vm(i) = vm(i)*(1 + c(i)*(ratio - 0.6)/0.4);
      elseif ratio < 0.4
         vm(i) = vm(i)/(1 + c(i)*((0.4 - ratio)/0.4));
      end
      if vm(i) > (ub(i)-lb(i))
         vm(i) = ub(i) - lb(i);
      end
   end

   %The pct. of times we did not find a move along a dimension 
   for i = 1:n
      pct_no_move(i) = round(num_no_move(i)/(ns*m)*100);
   end

   if iprint >= 2 
     prt8(n, vm, xopt, fopt, max, x, m);
     prtvec(pct_no_move, n, 'Pct of no moves');
     disp(['NUMBER OF FUNCTION EVALUATIONS ', num2str(nfcnev)]);
   end
  
   nacp(1:n) = 0;

   % Saving the results
   if isempty(Save_To_File) == 0
      save(Save_To_File);
   end

end %for m - the number of random walkers


if iprint >= 1
   prt9(max,n,t,xopt,vm,fopt,nup,ndown,nrej,lnobds,nnew)
   prtvec(pct_no_move, n, 'Pct of no moves')
   disp(['NUMBER OF FUNCTION EVALUATIONS  ', num2str(nfcnev)]);
end

%Check termination criteria.
quit = 0;
fstar(1) = f;
if (fopt - fstar(1)) <= eps 
    quit = 1;
end

for i=1:neps;
   if abs(f - fstar(i)) > eps
       quit = 0;
   end
end

%Terminate SA if appropriate.
if quit == 1
  x(1:n) = xopt(1:n);
  ier = 0;
  if max == 1
      fopt = fopt;
  else
      fopt = -fopt;
  end
  if iprint >= 0 
      %prt10()
      disp('***   SIMULATED ANNEALING HAS ACHIEVED SUCCESSFUL TERMINATION CRITERIA. IER = 0. ***');
  end
end

%If termination criteria is not met, prepare for another loop.
t = rt*t;
for i=neps:-1:2;
  fstar(i) = fstar(i-1);
end

f = fopt;
x(1:n) = xopt(1:n);

% Saving the results
if isempty(Save_To_File) == 0
    save(Save_To_File);
end

end % for the while loop index by quit'




% function which are used in the algorithm

function prt3(max, n, xp, x, f);
disp(' ');
prtvec(x, n, 'CURRENT X');
if max == 1 
    disp(['CURRENT F:' num2str(f)]);
else
    disp(['CURRENT F:' num2str(-f)]);
end
prtvec(xp, n, 'TRIAL X')
disp('POINT REJECTED SINCE OUT OF BOUNDS');


function prt4(max, n, xp, x, fp, f);
disp(' ');
prtvec(x,n,'CURRENT X');
if max == 1
  disp(['CURRENT F: ', num2str(f)]);
  prtvec(xp,n,'TRIAL X');
  disp(['RESULTING F: ', num2str(fp)]);
else
  disp(['CURRENT F: ', num2str(-f)]);
  prtvec(xp,n,'TRIAL X');
  disp(['RESULTING F: ', num2str(-fp)]);
end


function prt6(max);
if max == 1
  disp('THOUGH LOWER, POINT ACCEPTED');
else
  disp('THOUGH HIGHER, POINT ACCEPTED');
end


function prt7(max);
if max == 1
  disp('LOWER POINT REJECTED');
else
  disp('HIGHER POINT REJECTED');
end


function prt8(n, vm, xopt, fopt, max, x, m);
disp(' ');
disp('INTERMEDIATE RESULTS AFTER STEP LENGTH ADJUSTMENT')
disp(['AFTER RANDOM WALKER NUMBER ', num2str(m)]);
disp(' ');
prtvec(vm, n, 'NEW STEP LENGTH (VM)');
disp(' ');
prtvec(xopt, n, 'CURRENT OPTIMAL X');
disp(' ');
prtvec(x, n, 'CURRENT X');
disp(' ');
if max == 1
    disp(['CURRENT OPTIMAL VALUE OF THE OBJECT FUNCTION ', num2str(fopt)]);
else
    disp(['CURRENT OPTIMAL VALUE OF THE OBJECT FUNCTION ', num2str(-fopt)]);
end
disp(' ');


function prt9(max, n, t, xopt, vm, fopt, nup, ndown, nrej, lnobds, nnew)
totmov = nup + ndown + nrej;
  disp(' ');
  disp(' ');
  disp(' ');  
  disp('INTERMEDIATE RESULTS BEFORE NEXT TEMPERATURE REDUCTION');
  disp(' ');
  disp(['CURRENT TEMPERATURE:            ', num2str(t)]);
if max == 1
  disp(['MAX FUNCTION VALUE SO FAR:      ', num2str(fopt)]);
  disp(['TOTAL MOVES:                    ', num2str(totmov)]);
  disp(['UPHILL:                         ', num2str(nup)]);
  disp(['ACCEPTED DOWNHILL:              ', num2str(ndown)]);
  disp(['REJECTED DOWNHILL:              ', num2str(nrej)]);
  disp(['OUT OF BOUNDS TRIALS:           ', num2str(lnobds)]);
  disp(['NEW MAXIMA THIS TEMPERATURE:    ', num2str(nnew)]);
  disp(' ');
else
  disp(['MIN FUNCTION VALUE SO FAR:      ', num2str(-fopt)]);
  disp(['TOTAL MOVES:                    ', num2str(totmov)]);
  disp(['DOWNHILL:                       ', num2str(nup)]);
  disp(['ACCEPTED UPHILL:                ', num2str(ndown)]);
  disp(['REJECTED UPHILL:                ', num2str(nrej)]);
  disp(['TRIALS OUT OF BOUNDS:           ', num2str(lnobds)]);
  disp(['NEW MINIMA THIS TEMPERATURE:    ', num2str(nnew)]);
  disp(' ');
end
prtvec(xopt, n, 'CURRENT OPTIMAL X');
prtvec(vm, n, 'STEP LENGTH (VM)');
disp(' ');


function prtvec(vector_input, ncols, name);
%This subroutine prints the double precision vector named VECTOR.
%Elements 1 thru NCOLS will be printed. NAME is a character variable
%that describes VECTOR. Note that if NAME is given in the call to
%PRTVEC, it must be enclosed in quotes. If there are more than 5
%elements in VECTOR, 5 elements will be printed on each line.

% The current format
format_old = get(0,'Format');
% Format for displaying
format('short');
disp(name);
vector(1,:) = vector_input(:);
if ncols > 5
  lines = floor(ncols/5);
  for i = 1:lines
    ll = 5*(i - 1);
    disp(vector(1+ll:5+ll));
  end
  if ll+5 < ncols
      disp(vector(ll+5+1:ncols));
  end
else
  disp(vector(1:ncols));
end

% Resetting the format
format(format_old);



