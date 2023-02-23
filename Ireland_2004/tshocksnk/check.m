% check.m: Checks the solution found by solv.m for the New Keynesian model
%            with markup and technology shocks by making sure that the
%            equilibrium conditions hold after each type of shock. If the
%            solution is accurate, this program will return a 8x4 martix
%            of zeros.
%
% THIS PROGRAM WAS WRITTEN FOR MATLAB BY
%
%   PETER N. IRELAND
%   BOSTON COLLEGE
%   DEPARTMENT OF ECONOMICS
%   140 COMMONWEALTH AVENUE
%   CHESTNUT HILL, MA 02467
%   irelandp@bc.edu
%   http://www2.bc.edu/~irelandp
%
%  FINANCIAL SUPPORT FROM THE NATIONAL SCIENCE FOUNDATION UNDER GRANT NO.
%    SES-0213461 IS GRATEFULLY ACKNOWLEDGED.
%
%  COPYRIGHT (c) 2003 BY PETER N. IRELAND.  REDISTRIBUTION IS PERMITTED FOR
%    EDUCATIONAL AND RESEARCH PURPOSES, SO LONG AS NO CHANGES ARE MADE. ALL
%    COPIES MUST BE PROVIDED FREE OF CHARGE AND MUST INCLUDE THIS COPYRIGHT
%    NOTICE.

% check equilibrium conditions after an IS shock

  evec = [ siga 0 0 0 ]';

  st = bigw*evec;
  stp = bigpi*st;
  stpp = bigpi*stp;
  
  ylag = st(1);
  rlag = st(2);
  plag = st(3);
  glag = st(4);
  xlag = st(5);
  at = st(6);
  et = st(7);
  zt = st(8);
  epsrt = st(9);

  yt = stp(1);
  rt = stp(2);
  pt = stp(3);
  gt = stp(4);
  xt = stp(5);
  atp = stp(6);
  etp = stp(7);
  ztp = stp(8);
  epsrtp = stp(9);
  
  ytp = stpp(1);
  rtp = stpp(2);
  ptp = stpp(3);
  gtp = stpp(4);
  xtp = stpp(5);

  eqconda = ones(8,1);
  
  eqconda(1) = atp - rhoa*at;
  
  eqconda(2) = xt - alphax*xlag - (1-alphax)*xtp ...
                  + rt - ptp - (1-omega)*(1-rhoa)*at;
  
  eqconda(3) = etp - rhoe*et;
  
  eqconda(4) = ztp;
  
  eqconda(5) = pt - beta*alphap*plag - beta*(1-alphap)*ptp - psi*xt + et;
             
  eqconda(6) = gtp - ytp + yt - ztp;
  
  eqconda(7) = xt - yt + omega*at;
  
  eqconda(8) = rtp - rhor*rt - rhop*ptp - rhog*gtp - rhox*xtp;
  
% check equilibrium conditions after a markup shock

  evec = [ 0 sige 0 0 ]';

  st = bigw*evec;
  stp = bigpi*st;
  stpp = bigpi*stp;
  
  ylag = st(1);
  rlag = st(2);
  plag = st(3);
  glag = st(4);
  xlag = st(5);
  at = st(6);
  et = st(7);
  zt = st(8);
  epsrt = st(9);

  yt = stp(1);
  rt = stp(2);
  pt = stp(3);
  gt = stp(4);
  xt = stp(5);
  atp = stp(6);
  etp = stp(7);
  ztp = stp(8);
  epsrtp = stp(9);
  
  ytp = stpp(1);
  rtp = stpp(2);
  ptp = stpp(3);
  gtp = stpp(4);
  xtp = stpp(5);

  eqconde = ones(8,1);
  
  eqconde(1) = atp - rhoa*at;
  
  eqconde(2) = xt - alphax*xlag - (1-alphax)*xtp ...
                  + rt - ptp - (1-omega)*(1-rhoa)*at;
  
  eqconde(3) = etp - rhoe*et;
  
  eqconde(4) = ztp;
  
  eqconde(5) = pt - beta*alphap*plag - beta*(1-alphap)*ptp - psi*xt + et;
             
  eqconde(6) = gtp - ytp + yt - ztp;
  
  eqconde(7) = xt - yt + omega*at;
  
  eqconde(8) = rtp - rhor*rt - rhop*ptp - rhog*gtp - rhox*xtp;

% check equilibrium conditions after a technology shock

  evec = [ 0 0 sigz 0 ]';

  st = bigw*evec;
  stp = bigpi*st;
  stpp = bigpi*stp;
  
  ylag = st(1);
  rlag = st(2);
  plag = st(3);
  glag = st(4);
  xlag = st(5);
  at = st(6);
  et = st(7);
  zt = st(8);
  epsrt = st(9);

  yt = stp(1);
  rt = stp(2);
  pt = stp(3);
  gt = stp(4);
  xt = stp(5); 
  atp = stp(6);
  etp = stp(7);
  ztp = stp(8);
  epsrtp = stp(9);
  
  ytp = stpp(1);
  rtp = stpp(2);
  ptp = stpp(3);
  gtp = stpp(4);
  xtp = stpp(5);

  eqcondz = ones(8,1);
  
  eqcondz(1) = atp - rhoa*at;
  
  eqcondz(2) = xt - alphax*xlag - (1-alphax)*xtp ...
                  + rt - ptp - (1-omega)*(1-rhoa)*at;
  
  eqcondz(3) = etp - rhoe*et;
  
  eqcondz(4) = ztp;
  
  eqcondz(5) = pt - beta*alphap*plag - beta*(1-alphap)*ptp - psi*xt + et;
             
  eqcondz(6) = gtp - ytp + yt - ztp;
  
  eqcondz(7) = xt - yt + omega*at;
  
  eqcondz(8) = rtp - rhor*rt - rhop*ptp - rhog*gtp - rhox*xtp;
  
% check equilibrium conditions after a monetary policy shock

  evec = [ 0 0 0 sigr ]';

  st = bigw*evec;
  stp = bigpi*st;
  stpp = bigpi*stp;
  
  ylag = st(1);
  rlag = st(2);
  plag = st(3);
  glag = st(4);
  xlag = st(5);
  at = st(6);
  et = st(7);
  zt = st(8);
  epsrt = st(9);

  yt = stp(1);
  rt = stp(2);
  pt = stp(3);
  gt = stp(4);
  xt = stp(5);
  atp = stp(6);
  etp = stp(7);
  ztp = stp(8);
  epsrtp = stp(9);
  
  ytp = stpp(1);
  rtp = stpp(2);
  ptp = stpp(3);
  gtp = stpp(4);
  xtp = stpp(5);

  eqcondr = ones(8,1);
  
  eqcondr(1) = atp - rhoa*at;
  
  eqcondr(2) = xt - alphax*xlag - (1-alphax)*xtp ...
                  + rt - ptp - (1-omega)*(1-rhoa)*at;
  
  eqcondr(3) = etp - rhoe*et;
  
  eqcondr(4) = ztp;
  
  eqcondr(5) = pt - beta*alphap*plag - beta*(1-alphap)*ptp - psi*xt + et;
             
  eqcondr(6) = gtp - ytp + yt - ztp;
  
  eqcondr(7) = xt - yt + omega*at;
  
  eqcondr(8) = rtp - rhor*rt - rhop*ptp - rhog*gtp - rhox*xtp;
  
% report results

  [ eqconda eqconde eqcondz eqcondr ]