% stabtest.m: Performs parameter stability tests for the New Keynesian
%               model with markup and technology shocks.  Breakpoint occurs
%               at 1980:1.
%
%             Returns:
%
%               tstar1 = estimated parameters, 1948:1-1979:4
%               sevec1 = standard errors, 1948:1-1979:4
%               tstar2 = estimated parameters, 1980:1-2003:1
%               sevec2 = standard errors, 1980:1-2003:1
%               lrstat = LR statistic for stability of all parameters
%               wstat = Wald statistic for stability of all parameters
%               wstatp = Wald statistic for stability of policy parameters
%               wstatn = Wald statistic for stability of non-policy
%                          parameters
%               wstats = Wald statistic for stability of the standard
%                          deviations of the shocks
%               wstatvec = vector of Wald statistics for stability of each
%                            parameter
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

% load data and set pre-1980 subsample

  global gt pt rt scalinv

  load gpr.dat;
  
  gt = gpr(1:127,1);
  pt = gpr(1:127,2);
  rt = gpr(1:127,3);

% set starting values for pre-1980 subsample

  ztr = 0.0048;
  ptr = 0.0086;
  betatr = sqrt(0.99/0.01)/100;
  
  omegatr = 0.0000;
  psitr = 0.10;
  alphaxtr = 0.5000;
  alphaptr = 0.0000;

  rhortr = 1;
  rhoptr = 0.3000;
  rhogtr = 0.2500;
  rhoxtr = 0.0000;
  
  rhoatr = 0.1000;
  rhoetr = 0.0100;
  
  sigatr = 0.1500;
  sigetr = 0.0035;
  sigztr = 0.0100;
  sigrtr = 0.0030;
  
  bigtheto = [ omegatr alphaxtr alphaptr ...
               rhoptr rhogtr rhoxtr ...
               rhoatr rhoetr ...
               sigatr sigetr sigztr sigrtr ]';
              
  scalvec = [ 1 ; 1 ; 1 ; 1 ; 1; 100 ; 1 ; 100 ; ones(4,1) ];

% maximize likelihood for pre-1980 subsample

  options = optimset('Display','iter','LargeScale','off','MaxFunEvals',2500,'MaxIter',2500);

  thetstar = fminunc(@llfn,bigtheto,options);

% find standard errors for pre-1980 subsample

  thetstar = real(thetstar);

  omegatr = thetstar(1);
  alphaxtr = thetstar(2);
  alphaptr = thetstar(3);
  
  rhoptr = thetstar(4);
  rhogtr = thetstar(5);
  rhoxtr = thetstar(6);
  
  rhoatr = thetstar(7);
  rhoetr = thetstar(8);
  
  sigatr = thetstar(9);
  sigetr = thetstar(10);
  sigztr = thetstar(11);
  sigrtr = thetstar(12);
  
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
  
  tstar = [ omega alphax alphap ...
            rhop rhog rhox ...
            rhoa rhoe ...
            siga sige sigz sigr ]';

  scalinv = inv(diag(scalvec));

  tstars = diag(scalvec)*tstar;

  fstar = llfnse(tstars);

  eee = 1e-6;

  epsmat = eee*eye(12);

  hessvec = zeros(12,1);

  for i = 1:12

    hessvec(i) = llfnse(tstars+epsmat(:,i));

  end

  hessmatl = zeros(12,12);

  for i = 1:12

    for j = 1:i

      hessmatl(i,j) = (llfnse(tstars+epsmat(:,i)+epsmat(:,j)) ...
                        -hessvec(i)-hessvec(j)+fstar)/(eee^2);

    end

  end

  hessmatu = hessmatl' - diag(diag(hessmatl));

  hessmatf = hessmatl + hessmatu;

  bighx = scalinv*inv(hessmatf)*scalinv';

  sevec = sqrt(diag(bighx));

  tstar1 = tstar;

  bighx1 = bighx;

  sevec1 = sevec;
  
  fstar1 = fstar;

% set post-1980 subsample

  global gt pt rt scalinv

  load gpr.dat;
  
  gt = gpr(128:220,1);
  pt = gpr(128:220,2);
  rt = gpr(128:220,3);
  
% set starting values for post-1980 subsample

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
 
  scalvec = [ 1 ; 100 ; 100 ; ones(9,1) ];

% maximize likelihood for post-1980 subsample

  options = optimset('Display','iter','LargeScale','off','MaxFunEvals',2500,'MaxIter',2500);

  thetstar = fminunc(@llfn,bigtheto,options);

% find standard errors for post-1980 subsample

  thetstar = real(thetstar);

  omegatr = thetstar(1);
  alphaxtr = thetstar(2);
  alphaptr = thetstar(3);
  
  rhoptr = thetstar(4);
  rhogtr = thetstar(5);
  rhoxtr = thetstar(6);
  
  rhoatr = thetstar(7);
  rhoetr = thetstar(8);
  
  sigatr = thetstar(9);
  sigetr = thetstar(10);
  sigztr = thetstar(11);
  sigrtr = thetstar(12);
  
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
  
  tstar = [ omega alphax alphap ...
            rhop rhog rhox ...
            rhoa rhoe ...
            siga sige sigz sigr ]';

  scalinv = inv(diag(scalvec));

  tstars = diag(scalvec)*tstar;

  fstar = llfnse(tstars);

  eee = 1e-6;

  epsmat = eee*eye(12);

  hessvec = zeros(12,1);

  for i = 1:12

    hessvec(i) = llfnse(tstars+epsmat(:,i));

  end

  hessmatl = zeros(12,12);

  for i = 1:12

    for j = 1:i

      hessmatl(i,j) = (llfnse(tstars+epsmat(:,i)+epsmat(:,j)) ...
                        -hessvec(i)-hessvec(j)+fstar)/(eee^2);

    end

  end

  hessmatu = hessmatl' - diag(diag(hessmatl));

  hessmatf = hessmatl + hessmatu;

  bighx = scalinv*inv(hessmatf)*scalinv';

  sevec = sqrt(diag(bighx));

  tstar2 = tstar;

  bighx2 = bighx;

  sevec2 = sevec;
  
  fstar2 = fstar;

% form test statistics

  lrstat = -2*(fstar1+fstar2+2648.4303);

  wstat = (tstar1-tstar2)'*inv(bighx1+bighx2)*(tstar1-tstar2);

  bighxp1 = bighx1(4:6,4:6);
      
  bighxp2 = bighx2(4:6,4:6);
      
  tstarp1 = tstar1(4:6);
  
  tstarp2 = tstar2(4:6);

  wstatp = (tstarp1-tstarp2)'*inv(bighxp1+bighxp2)*(tstarp1-tstarp2);
        
  bighxn1 = [ bighx1(1:3,1:3) bighx1(1:3,7:8) ; ...
              bighx1(7:8,1:3) bighx1(7:8,7:8) ];
      
  bighxn2 = [ bighx2(1:3,1:3) bighx2(1:3,7:8) ; ...
              bighx2(7:8,1:3) bighx2(7:8,7:8) ];    
  
  tstarn1 = [ tstar1(1:3) ; tstar1(7:8) ];
  
  tstarn2 = [ tstar2(1:3) ; tstar2(7:8) ];

  wstatn = (tstarn1-tstarn2)'*inv(bighxn1+bighxn2)*(tstarn1-tstarn2);
  
  bighxs1 = bighx1(9:12,9:12);
  
  bighxs2 = bighx2(9:12,9:12);
  
  tstars1 = tstar1(9:12);
  
  tstars2 = tstar2(9:12);
  
  wstats = (tstars1-tstars2)'*inv(bighxs1+bighxs2)*(tstars1-tstars2);
  
  wstatvec = (tstar1-tstar2).*(tstar1-tstar2)./(diag(bighx1)+diag(bighx2));