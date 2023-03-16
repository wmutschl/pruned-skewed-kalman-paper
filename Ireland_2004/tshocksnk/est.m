% est.m: Maximizes the (minimizes the negative) log likelihood function for
%          the New Keynesian model with markup and technology shocks.
%
%        When maximizing the log likelihood function, the parameters are
%          transformed to satisfy theoretical restrictions. The log
%          likelihood function with transformed parameters is contained in
%          llfn.m.
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

% load data and choose sample

  global gt pt rt

  load gpr.dat;

%  gt = gpr(:,1);
%  pt = gpr(:,2);
%  rt = gpr(:,3);

%  gt = gpr(1:127,1);
%  pt = gpr(1:127,2);
%  rt = gpr(1:127,3);

  gt = gpr(128:220,1);
  pt = gpr(128:220,2);
  rt = gpr(128:220,3);
  
% set starting values for full sample

%  ztr = 0.0048;
%  ptr = 0.0086;
%  betatr = sqrt(0.99/0.01)/100;
  
%  omegatr = sqrt(0.625/0.375);
%  psitr = 0.10;
%  alphaxtr = sqrt(0.25/0.75);
%  alphaptr = sqrt(0.25/0.75);
  
%  rhortr = 1;
%  rhoptr = 0.3515;
%  rhogtr = 0.2556;
%  rhoxtr = 0.0285;
  
%  rhoatr = 0.0287;
%  rhoetr = 0.0676;
  
%  sigatr = 0.0187;
%  sigetr = 0.0014;
%  sigztr = 0.0067;
%  sigrtr = 0.0032;
  
%  bigtheto = [ omegatr alphaxtr alphaptr ...
%               rhoptr rhogtr rhoxtr ...
%               rhoatr rhoetr ...
%               sigatr sigetr sigztr sigrtr ]';
  
% set starting values for pre-1980 sample

%  ztr = 0.0048;
%  ptr = 0.0086;
%  betatr = sqrt(0.99/0.01)/100;
  
%  omegatr = 0.0000;
%  psitr = 0.10;
%  alphaxtr = 0.5000;
%  alphaptr = 0.0000;

%  rhortr = 1;
%  rhoptr = 0.3000;
%  rhogtr = 0.2500;
%  rhoxtr = 0.0000;
  
%  rhoatr = 0.1000;
%  rhoetr = 0.0100;
  
%  sigatr = 0.1500;
%  sigetr = 0.0035;
%  sigztr = 0.0100;
%  sigrtr = 0.0030;
  
%  bigtheto = [ omegatr alphaxtr alphaptr ...
%               rhoptr rhogtr rhoxtr ...
%               rhoatr rhoetr ...
%               sigatr sigetr sigztr sigrtr ]';
  
% set starting values for post-1980 sample

  ztr = 0.0048;
  ptr = 0.0086;
  betatr = sqrt(0.99/0.01)/100;
 
  omegatr = sqrt(0.625/0.375);
  psitr = 0.10;
  alphaxtr = sqrt(0.25/0.75);
  alphaptr = sqrt(0.25/0.75);
  
  rhortr = 1; 
  rhoptr = 0.6017;
  rhogtr = 0.4240;
  rhoxtr = 0.0770;
  
  rhoatr = 0.0205;
  rhoetr = 0.0882;
 
  sigatr = 0.0150;
  sigetr = 0.0008;
  sigztr = 0.0000;
  sigrtr = 0.0030;
  
  bigtheto = [ omegatr alphaxtr alphaptr ...
               rhoptr rhogtr rhoxtr ...
               rhoatr rhoetr ...
               sigatr sigetr sigztr sigrtr ]';

% maximize likelihood

  options = optimset('Display','iter','LargeScale','off','MaxFunEvals',10000,'MaxIter',10000);

  [thetstar,fval] = fminsearch(@llfn,bigtheto,options);

  omegatr = thetstar(1);
  omega = omegatr^2/(1+omegatr^2);
  alphaxtr = thetstar(2);
  alphax = alphaxtr^2/(1+alphaxtr^2);
  alphaptr = thetstar(3);
  alphap = alphaptr^2/(1+alphaptr^2);
  rhoptr = thetstar(4);
  rhop = abs(rhoptr);
  rhogtr = thetstar(5);
  rhog = abs(rhogtr);
  rhoxtr = thetstar(6);
  rhox = abs(rhoxtr);
  rhoatr = thetstar(7);
  rhoa = (100*rhoatr)^2/(1+(100*rhoatr)^2);
  rhoetr = thetstar(8);
  rhoe = (100*rhoetr)^2/(1+(100*rhoetr)^2);
  sigatr = thetstar(9);
  siga = abs(sigatr);
  sigetr = thetstar(10);
  sige = abs(sigetr);
  sigztr = thetstar(11);
  sigz = abs(sigztr);
  sigrtr = thetstar(12);
  sigr = abs(sigrtr);
  bigthetopt = [ omega alphax alphap ...
               rhop rhog rhox ...
               rhoa rhoe ...
               siga sige sigz sigr ]';
  fval1 = llfnse(bigthetopt);