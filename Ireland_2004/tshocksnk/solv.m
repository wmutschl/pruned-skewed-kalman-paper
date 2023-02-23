% solv.m: Solves the New Keynesian model with markup and technology shocks.
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

% set parameter values equal to full sample estimates

%  z = 1.0048;
%  p = 1.0086;
%  beta = 0.99;

%  omega = 0.0617;
%  psi = 0.10;
%  alphax = 0.0836;
%  alphap = 0;
  
%  rhor = 1;
%  rhop = 0.3597;
%  rhog = 0.2536;
%  rhox = 0.0347;
  
%  rhoa = 0.9470;
%  rhoe = 0.9625;

%  siga = 0.0405;
%  sige = 0.0012;
%  sigz = 0.0109;
%  sigr = 0.0031;

% set parameter values equal to pre-1980 estimates

%  z = 1.0048;
%  p = 1.0086;
%  beta = 0.99;

%  omega = 0;
%  psi = 0.10;
%  alphax = 0.2028;
%  alphap = 0;
  
%  rhor = 1;
%  rhop = 0.3053;
%  rhog = 0.2365;
%  rhox = 0;
  
%  rhoa = 0.9910;
%  rhoe = 0.5439;

%  siga = 0.1538;
%  sige = 0.0035;
%  sigz = 0.0104;
%  sigr = 0.0033;

% set parameter values equal to post-1980 estimates

  z = 1.0048;
  p = 1.0086;
  beta = 0.99;

  omega = 0.0581;
  psi = 0.10;
  alphax = 0;
  alphap = 0;
  
  rhor = 1;
  rhop = 0.3866;
  rhog = 0.3960;
  rhox = 0.1654;
  
  rhoa = 0.9048;
  rhoe = 0.9907;

  siga = 0.0302;
  sige = 0.0002;
  sigz = 0.0089;
  sigr = 0.0028;
  
% find steady state values

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

  if abs(bigt11(5,5)/bigs11(5,5)) > 1
    'error - no solution'
  end

  if abs(bigt22(1,1)/bigs22(1,1)) < 1
    'error - multiple solutions'
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