% imp.m: Uses the solution found by solv.m to compute impulse responses for
%          the New Keynesian model with markup and technology shocks.
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

% IS shock

  sa = zeros(51,9);

  evec = [ 100*siga 0 0 0 ]';

  sa(2,:) = (bigw*evec)';

  for t = 3:51

    sa(t,:) = (bigpi*sa(t-1,:)')';

  end

% markup shock

  se = zeros(51,9);

  evec = [ 0 100*sige 0 0 ]';

  se(2,:) = (bigw*evec)';

  for t = 3:51

    se(t,:) = (bigpi*se(t-1,:)')';

  end

% technology shock

  sz = zeros(51,9);

  evec = [ 0 0 100*sigz 0 ]';

  sz(2,:) = (bigw*evec)';
  
  for t = 3:51

    sz(t,:) = (bigpi*sz(t-1,:)')';

  end

% monetary policy shock

  sr = zeros(51,9);
  
  evec = [ 0 0 0 100*sigr ]';

  sr(2,:) = (bigw*evec)';

  for t = 3:51

    sr(t,:) = (bigpi*sr(t-1,:)')';

  end

% create output vectors

  ya = sa(2:51,1);
  ra = sa(2:51,2);
  pa = sa(2:51,3);
  ga = sa(2:51,4);
  xa = sa(2:51,5);
  aa = sa(1:50,6);
  ea = sa(1:50,7);
  za = sa(1:50,8);
  epsra = sa(1:50,9);
  
  ye = se(2:51,1);
  re = se(2:51,2);
  pe = se(2:51,3);
  ge = se(2:51,4);
  xe = se(2:51,5);
  ae = se(1:50,6);
  ee = se(1:50,7);
  ze = se(1:50,8);
  epsre = se(1:50,9);
  
  yz = sz(2:51,1);
  rz = sz(2:51,2);
  pz = sz(2:51,3);
  gz = sz(2:51,4);
  xz = sz(2:51,5);
  az = sz(1:50,6);
  ez = sz(1:50,7);
  zz = sz(1:50,8);
  epsrz = sz(1:50,9);
  
  yr = sr(2:51,1);
  rr = sr(2:51,2);
  pr = sr(2:51,3);
  gr = sr(2:51,4);
  xr = sr(2:51,5);
  ar = sr(1:50,6);
  er = sr(1:50,7);
  zr = sr(1:50,8);
  epsrr = sr(1:50,9);