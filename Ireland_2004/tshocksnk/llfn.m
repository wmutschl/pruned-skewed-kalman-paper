function llfn = llfn(bigthet);
% Uses the Kalman filter to evaluate the negative log likelihood function
%   for the New Keynesian model with markup and technology shocks. The
%   parameters are transformed to satisfy theoretical restrictions.
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

  global gt pt rt

  capt = length(gt);

  bigthet = real(bigthet);
  
  ztr = 0.0048;
  ptr = 0.0086;
  betatr = sqrt(0.99/0.01)/100;

  omegatr = sqrt(0.625/0.375);
  psitr = 0.10;
  alphaxtr = sqrt(0.25/0.75);
  alphaptr = sqrt(0.25/0.75);
  
  rhortr = 1;
  rhoptr = 0.25;
  rhogtr = 0.05;
  rhoxtr = 0.05;
  
  rhoatr = 0;
  rhoetr = 0;

  sigatr = 0.1;
  sigetr = 0.1;
  sigztr = 0.1;
  sigrtr = 0.1;
  
  omegatr = bigthet(1);
  
  alphaxtr = bigthet(2);
  alphaptr = bigthet(3);
  
  rhoptr = bigthet(4);
  rhogtr = bigthet(5);
  rhoxtr = bigthet(6);
  
  rhoatr = bigthet(7);
  rhoetr = bigthet(8);
  
  sigatr = bigthet(9);
  sigetr = bigthet(10);
  sigztr = bigthet(11);
  sigrtr = bigthet(12);

% untransform parameters

  z = 1 + abs(ztr);
  p = 1 + abs(ptr);
  beta = (100*betatr)^2/(1+(100*betatr)^2);

  omega = omegatr^2/(1+omegatr^2);
  psi = abs(psitr);
  alphax = alphaxtr^2/(1+alphaxtr^2);
  alphap = alphaptr^2/(1+alphaptr^2);
  
  rhor = abs(rhortr);
  rhop = abs(rhoptr);
  rhog = abs(rhogtr);
  rhox = abs(rhoxtr);
  
  rhoa = (100*rhoatr)^2/(1+(100*rhoatr)^2);
  rhoe = (100*rhoetr)^2/(1+(100*rhoetr)^2);
  
  siga = abs(sigatr);
  sige = abs(sigetr);
  sigz = abs(sigztr);
  sigr = abs(sigrtr);
  
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

  llfn = (3*capt/2)*log(2*pi);

  for t = 1:capt

    ut = dthat(t,:)' - bigcx*st;

    omegt = bigcx*bigsigt*bigcx';

    omeginvt = inv(omegt);

    llfn = llfn + (1/2)*(log(det(omegt))+ut'*omeginvt*ut);

    bigkt = bigax*bigsigt*bigcx'*omeginvt;

    st = bigax*st + bigkt*ut;

    bigsigt = bigbvbx + bigax*bigsigt*bigax' ...
                - bigax*bigsigt*bigcx'*omeginvt*bigcx*bigsigt*bigax';

  end

% penalize eigenvalue constraint violations

  if lamviol

    llfn = llfn + 1e12;

  end

  if abs(imag(llfn)) > 0

    llfn = real(llfn) + 1e12;

  end