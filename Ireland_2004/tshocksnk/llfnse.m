function llfnse = llfnse(bigthet);
% Uses the Kalman filter to evaluate the negative log likelihood function
%   for the New Keynesian model with markup and technology shocks. The
%   parameters are untransformed to facilitate calculation of standard
%   errors.
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

% define variables and parameters

  global gt pt rt scalinv

  bigthet = scalinv*bigthet;

  capt = length(gt);
  
  z = 1.0048;
  p = 1.0086;
  beta = 0.99;

  omega = bigthet(1);
  psi = 0.10;
  alphax = bigthet(2);
  alphap = bigthet(3);
  
  rhor = 1;
  rhop = bigthet(4);
  rhog = bigthet(5);
  rhox = bigthet(6);
  
  rhoa = bigthet(7);
  rhoe = bigthet(8);

  siga = bigthet(9);
  sige = bigthet(10);
  sigz = bigthet(11);
  sigr = bigthet(12);
  
% find steady state

  gss = z;  

  pss = p;
  
  rss = pss*(z/beta);

% form matrices A, B, and C

  biga = [ 0 -1 0 0 0 1 1-alphax ; ...
           0 0 0 0 psi beta*(1-alphap) 0 ; ...
           -1 0 0 1 0 0 0 ; ...
           1 0 0 0 0 0 0 ; ...
           0 1 -rhop -rhog -rhox 0 0 ; ...
           0 0 1 0 0 0 0 ; ...
           0 0 0 0 1 0 0 ];
   
   bigb = [ 0 0 0 0 -alphax 0 1 ; ...
            0 0 -beta*alphap 0 0 1 0 ; ...
            -1 0 0 0 0 0 0 ; ...
            0 0 0 0 0 0 1 ; ...
            0 rhor 0 0 0 0 0 ; ...
            0 0 0 0 0 1 0 ; ...
            0 0 0 0 0 0 1 ];
    
    bigc = [ (omega-1)*(1-rhoa) 0 0 0 ; ...
             0 1 0 0 ; ...
             0 0 1 0 ; ...
             omega 0 0 0 ; ...
             0 0 0 1 ; ...
             0 0 0 0 ;
             0 0 0 0 ];
  
% form matrix P

  bigp = [ rhoa 0 0 0 ; ...
           0 rhoe 0 0 ; ...
           0 0 0 0 ; ...
           0 0 0 0 ];
  
% form matrices Q, Z, S, and T  
  
  [bigs,bigt,bigq,bigz] = qz(biga,bigb);
  
  [bigs,bigt,bigq,bigz] = reorder(bigs,bigt,bigq,bigz);

  bigq1 = bigq(1:5,:);
  bigq2 = bigq(6:7,:);
  
  bigz11 = bigz(1:5,1:5);
  bigz12 = bigz(1:5,6:7);
  bigz21 = bigz(6:7,1:5);
  bigz22 = bigz(6:7,6:7);
  
  bigs11 = bigs(1:5,1:5);
  bigs12 = bigs(1:5,6:7);
  bigs22 = bigs(6:7,6:7);

  bigt11 = bigt(1:5,1:5);
  bigt12 = bigt(1:5,6:7);
  bigt22 = bigt(6:7,6:7);
  
  lamviol = 0;

  if abs(bigt11(5,5)/bigs11(5,5)) > 1
    lamviol = 1;
  end

  if abs(bigt22(1,1)/bigs22(1,1)) < 1
    lamviol = 1;
  end

% form matrix R

  bigra = bigs22*inv(bigt22);
  
  bigrb = bigq2*bigc;
  
  vecr = inv(eye(8)-kron(bigp,bigra))*bigrb(:);

  bigr = reshape(vecr,2,4);

% form matrices M

  bigm3 = bigz11*inv(bigs11)*bigt11*inv(bigz11);
  
  bigm4a = bigt11*inv(bigz11)*bigz12*inv(bigt22)*bigr ...
             + bigq1*bigc + bigs12*inv(bigt22)*bigr*bigp ...
             - bigt12*inv(bigt22)*bigr;
  
  bigm4 = bigz11*inv(bigs11)*bigm4a - bigz12*inv(bigt22)*bigr*bigp;

% form matrices PI and W

  bigpi = [ bigm3 bigm4 ; zeros(4,5) bigp ];

  bigw = [ zeros(5,4) ; eye(4) ];

  bigpi = real(bigpi);
    
% form matrices AX, BX, CX, VX, and BVBX

  bigax = bigpi;

  bigbx = bigw;

  bigcx = [ bigpi(4,:) ; bigpi(3,:) ; bigpi(2,:) ];

  bigvx = diag([siga^2 sige^2 sigz^2 sigr^2]);

  bigbvbx = bigbx*bigvx*bigbx';

% put data in deviation form

%  gthat = gt - log(gss);

%  pthat = pt - log(pss);

%  rthat = rt - log(rss);

  gthat = gt - mean(gt);
  
  pthat = pt - mean(pt);
  
  rthat = rt - mean(rt);

  dthat = [ gthat pthat rthat ];

% evaluate negative log likelihood

  st = zeros(9,1);

  bigsig1 = inv(eye(81)-kron(bigax,bigax))*bigbvbx(:);

  bigsigt = reshape(bigsig1,9,9);

  llfnse = (3*capt/2)*log(2*pi);

  for t = 1:capt

    ut = dthat(t,:)' - bigcx*st;

    omegt = bigcx*bigsigt*bigcx';

    omeginvt = inv(omegt);

    llfnse = llfnse + (1/2)*(log(det(omegt))+ut'*omeginvt*ut);

    bigkt = bigax*bigsigt*bigcx'*omeginvt;

    st = bigax*st + bigkt*ut;

    bigsigt = bigbvbx + bigax*bigsigt*bigax' ...
                - bigax*bigsigt*bigcx'*omeginvt*bigcx*bigsigt*bigax';

  end

% penalize eigenvalue constraint violations

  if lamviol

    llfnse = llfnse + 1e12;

  end

  if abs(imag(llfnse)) > 0

    llfnse = real(llfnse) + 1e12;

  end