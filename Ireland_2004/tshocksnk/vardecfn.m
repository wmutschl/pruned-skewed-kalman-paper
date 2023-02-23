function vardecfn = vardecfn(bigthet);
% For given parameter values, computes variance decompositions for the New
%   Keynesian model with markup and technology shocks.
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

  bigcx = [ bigpi(4,:) ; bigpi(3,:) ; bigpi(2,:) ; bigpi(5,:) ];

  bigvx = diag([siga^2 sige^2 sigz^2 sigr^2]);

  bigbvbx = bigbx*bigvx*bigbx';

% calculate k-step ahead forecast error variances

  vardecfn = zeros(28,5);

  bigfx = eye(9);

  bigsigsk = zeros(9,9);

  for k = 1:40

    bigsigsk = bigsigsk + bigfx*bigbvbx*bigfx';

    bigsigdk = bigcx*bigsigsk*bigcx';

    if k == 1
      vardecfn(1,1) = bigsigdk(1,1);
      vardecfn(8,1) = bigsigdk(2,2);
      vardecfn(15,1) = bigsigdk(3,3);
      vardecfn(22,1) = bigsigdk(4,4);
    end

    if k == 4
      vardecfn(2,1) = bigsigdk(1,1);
      vardecfn(9,1) = bigsigdk(2,2);
      vardecfn(16,1) = bigsigdk(3,3);
      vardecfn(23,1) = bigsigdk(4,4);
    end

    if k == 8
      vardecfn(3,1) = bigsigdk(1,1);
      vardecfn(10,1) = bigsigdk(2,2);
      vardecfn(17,1) = bigsigdk(3,3);
      vardecfn(24,1) = bigsigdk(4,4);
    end

    if k == 12
      vardecfn(4,1) = bigsigdk(1,1);
      vardecfn(11,1) = bigsigdk(2,2);
      vardecfn(18,1) = bigsigdk(3,3);
      vardecfn(25,1) = bigsigdk(4,4);
    end

    if k == 20
      vardecfn(5,1) = bigsigdk(1,1);
      vardecfn(12,1) = bigsigdk(2,2);
      vardecfn(19,1) = bigsigdk(3,3);
      vardecfn(26,1) = bigsigdk(4,4);
    end

    if k == 40
      vardecfn(6,1) = bigsigdk(1,1);
      vardecfn(13,1) = bigsigdk(2,2);
      vardecfn(20,1) = bigsigdk(3,3);
      vardecfn(27,1) = bigsigdk(4,4);
    end

    bigfx = bigfx*bigax;

  end    

  bigsigs1 = inv(eye(81)-kron(bigax,bigax))*bigbvbx(:);

  bigsigs = reshape(bigsigs1,9,9);

  bigsigd = bigcx*bigsigs*bigcx';

  vardecfn(7,1) = bigsigd(1,1);
  vardecfn(14,1) = bigsigd(2,2);
  vardecfn(21,1) = bigsigd(3,3);
  vardecfn(28,1) = bigsigdk(4,4);

% calculate percentages due to IS shocks

  bigvx = diag([siga^2 0 0 0]);

  bigbvbx = bigbx*bigvx*bigbx';

  bigfx = eye(9);

  bigsigsk = zeros(9,9);

  for k = 1:40

    bigsigsk = bigsigsk + bigfx*bigbvbx*bigfx';

    bigsigdk = bigcx*bigsigsk*bigcx';

    if k == 1
      vardecfn(1,2) = bigsigdk(1,1)/vardecfn(1,1);
      vardecfn(8,2) = bigsigdk(2,2)/vardecfn(8,1);
      vardecfn(15,2) = bigsigdk(3,3)/vardecfn(15,1);
      vardecfn(22,2) = bigsigdk(4,4)/vardecfn(22,1);
    end

    if k == 4
      vardecfn(2,2) = bigsigdk(1,1)/vardecfn(2,1);
      vardecfn(9,2) = bigsigdk(2,2)/vardecfn(9,1);
      vardecfn(16,2) = bigsigdk(3,3)/vardecfn(16,1);
      vardecfn(23,2) = bigsigdk(4,4)/vardecfn(23,1);
    end

    if k == 8
      vardecfn(3,2) = bigsigdk(1,1)/vardecfn(3,1);
      vardecfn(10,2) = bigsigdk(2,2)/vardecfn(10,1);
      vardecfn(17,2) = bigsigdk(3,3)/vardecfn(17,1);
      vardecfn(24,2) = bigsigdk(4,4)/vardecfn(24,1);
    end

    if k == 12
      vardecfn(4,2) = bigsigdk(1,1)/vardecfn(4,1);
      vardecfn(11,2) = bigsigdk(2,2)/vardecfn(11,1);
      vardecfn(18,2) = bigsigdk(3,3)/vardecfn(18,1);
      vardecfn(25,2) = bigsigdk(4,4)/vardecfn(25,1);
    end

    if k == 20
      vardecfn(5,2) = bigsigdk(1,1)/vardecfn(5,1);
      vardecfn(12,2) = bigsigdk(2,2)/vardecfn(12,1);
      vardecfn(19,2) = bigsigdk(3,3)/vardecfn(19,1);
      vardecfn(26,2) = bigsigdk(4,4)/vardecfn(26,1);
    end

    if k == 40
      vardecfn(6,2) = bigsigdk(1,1)/vardecfn(6,1);
      vardecfn(13,2) = bigsigdk(2,2)/vardecfn(13,1);
      vardecfn(20,2) = bigsigdk(3,3)/vardecfn(20,1);
      vardecfn(27,2) = bigsigdk(4,4)/vardecfn(27,1);
    end

    bigfx = bigfx*bigax;

  end    

  bigsigs1 = inv(eye(81)-kron(bigax,bigax))*bigbvbx(:);

  bigsigs = reshape(bigsigs1,9,9);

  bigsigd = bigcx*bigsigs*bigcx';

  vardecfn(7,2) = bigsigd(1,1)/vardecfn(7,1);
  vardecfn(14,2) = bigsigd(2,2)/vardecfn(14,1);
  vardecfn(21,2) = bigsigd(3,3)/vardecfn(21,1);
  vardecfn(28,2) = bigsigdk(4,4)/vardecfn(28,1);

% calculate percentages due to markup shocks

  bigvx = diag([0 sige^2 0 0]);

  bigbvbx = bigbx*bigvx*bigbx';

  bigfx = eye(9);

  bigsigsk = zeros(9,9);

  for k = 1:40

    bigsigsk = bigsigsk + bigfx*bigbvbx*bigfx';

    bigsigdk = bigcx*bigsigsk*bigcx';

    if k == 1
      vardecfn(1,3) = bigsigdk(1,1)/vardecfn(1,1);
      vardecfn(8,3) = bigsigdk(2,2)/vardecfn(8,1);
      vardecfn(15,3) = bigsigdk(3,3)/vardecfn(15,1);
      vardecfn(22,3) = bigsigdk(4,4)/vardecfn(22,1);
    end

    if k == 4
      vardecfn(2,3) = bigsigdk(1,1)/vardecfn(2,1);
      vardecfn(9,3) = bigsigdk(2,2)/vardecfn(9,1);
      vardecfn(16,3) = bigsigdk(3,3)/vardecfn(16,1);
      vardecfn(23,3) = bigsigdk(4,4)/vardecfn(23,1);
    end

    if k == 8
      vardecfn(3,3) = bigsigdk(1,1)/vardecfn(3,1);
      vardecfn(10,3) = bigsigdk(2,2)/vardecfn(10,1);
      vardecfn(17,3) = bigsigdk(3,3)/vardecfn(17,1);
      vardecfn(24,3) = bigsigdk(4,4)/vardecfn(24,1);
    end

    if k == 12
      vardecfn(4,3) = bigsigdk(1,1)/vardecfn(4,1);
      vardecfn(11,3) = bigsigdk(2,2)/vardecfn(11,1);
      vardecfn(18,3) = bigsigdk(3,3)/vardecfn(18,1);
      vardecfn(25,3) = bigsigdk(4,4)/vardecfn(25,1);
    end

    if k == 20
      vardecfn(5,3) = bigsigdk(1,1)/vardecfn(5,1);
      vardecfn(12,3) = bigsigdk(2,2)/vardecfn(12,1);
      vardecfn(19,3) = bigsigdk(3,3)/vardecfn(19,1);
      vardecfn(26,3) = bigsigdk(4,4)/vardecfn(26,1);
    end

    if k == 40
      vardecfn(6,3) = bigsigdk(1,1)/vardecfn(6,1);
      vardecfn(13,3) = bigsigdk(2,2)/vardecfn(13,1);
      vardecfn(20,3) = bigsigdk(3,3)/vardecfn(20,1);
      vardecfn(27,3) = bigsigdk(4,4)/vardecfn(27,1);
    end

    bigfx = bigfx*bigax;

  end    

  bigsigs1 = inv(eye(81)-kron(bigax,bigax))*bigbvbx(:);

  bigsigs = reshape(bigsigs1,9,9);

  bigsigd = bigcx*bigsigs*bigcx';

  vardecfn(7,3) = bigsigd(1,1)/vardecfn(7,1);
  vardecfn(14,3) = bigsigd(2,2)/vardecfn(14,1);
  vardecfn(21,3) = bigsigd(3,3)/vardecfn(21,1);
  vardecfn(28,3) = bigsigdk(4,4)/vardecfn(28,1);

% calculate percentages due to technology shocks

  bigvx = diag([0 0 sigz^2 0]);

  bigbvbx = bigbx*bigvx*bigbx';

  bigfx = eye(9);

  bigsigsk = zeros(9,9);

  for k = 1:40

    bigsigsk = bigsigsk + bigfx*bigbvbx*bigfx';

    bigsigdk = bigcx*bigsigsk*bigcx';

    if k == 1
      vardecfn(1,4) = bigsigdk(1,1)/vardecfn(1,1);
      vardecfn(8,4) = bigsigdk(2,2)/vardecfn(8,1);
      vardecfn(15,4) = bigsigdk(3,3)/vardecfn(15,1);
      vardecfn(22,4) = bigsigdk(4,4)/vardecfn(22,1);
    end

    if k == 4
      vardecfn(2,4) = bigsigdk(1,1)/vardecfn(2,1);
      vardecfn(9,4) = bigsigdk(2,2)/vardecfn(9,1);
      vardecfn(16,4) = bigsigdk(3,3)/vardecfn(16,1);
      vardecfn(23,4) = bigsigdk(4,4)/vardecfn(23,1);
    end

    if k == 8
      vardecfn(3,4) = bigsigdk(1,1)/vardecfn(3,1);
      vardecfn(10,4) = bigsigdk(2,2)/vardecfn(10,1);
      vardecfn(17,4) = bigsigdk(3,3)/vardecfn(17,1);
      vardecfn(24,4) = bigsigdk(4,4)/vardecfn(24,1);
    end

    if k == 12
      vardecfn(4,4) = bigsigdk(1,1)/vardecfn(4,1);
      vardecfn(11,4) = bigsigdk(2,2)/vardecfn(11,1);
      vardecfn(18,4) = bigsigdk(3,3)/vardecfn(18,1);
      vardecfn(25,4) = bigsigdk(4,4)/vardecfn(25,1);
    end

    if k == 20
      vardecfn(5,4) = bigsigdk(1,1)/vardecfn(5,1);
      vardecfn(12,4) = bigsigdk(2,2)/vardecfn(12,1);
      vardecfn(19,4) = bigsigdk(3,3)/vardecfn(19,1);
      vardecfn(26,4) = bigsigdk(4,4)/vardecfn(26,1);
    end

    if k == 40
      vardecfn(6,4) = bigsigdk(1,1)/vardecfn(6,1);
      vardecfn(13,4) = bigsigdk(2,2)/vardecfn(13,1);
      vardecfn(20,4) = bigsigdk(3,3)/vardecfn(20,1);
      vardecfn(27,4) = bigsigdk(4,4)/vardecfn(27,1);
    end

    bigfx = bigfx*bigax;

  end    

  bigsigs1 = inv(eye(81)-kron(bigax,bigax))*bigbvbx(:);

  bigsigs = reshape(bigsigs1,9,9);

  bigsigd = bigcx*bigsigs*bigcx';

  vardecfn(7,4) = bigsigd(1,1)/vardecfn(7,1);
  vardecfn(14,4) = bigsigd(2,2)/vardecfn(14,1);
  vardecfn(21,4) = bigsigd(3,3)/vardecfn(21,1);
  vardecfn(28,4) = bigsigdk(4,4)/vardecfn(28,1);

% calculate percentages due to policy shocks

  bigvx = diag([0 0 0 sigr^2]);

  bigbvbx = bigbx*bigvx*bigbx';

  bigfx = eye(9);

  bigsigsk = zeros(9,9);

  for k = 1:40

    bigsigsk = bigsigsk + bigfx*bigbvbx*bigfx';

    bigsigdk = bigcx*bigsigsk*bigcx';

    if k == 1
      vardecfn(1,5) = bigsigdk(1,1)/vardecfn(1,1);
      vardecfn(8,5) = bigsigdk(2,2)/vardecfn(8,1);
      vardecfn(15,5) = bigsigdk(3,3)/vardecfn(15,1);
      vardecfn(22,5) = bigsigdk(4,4)/vardecfn(22,1);
    end

    if k == 4
      vardecfn(2,5) = bigsigdk(1,1)/vardecfn(2,1);
      vardecfn(9,5) = bigsigdk(2,2)/vardecfn(9,1);
      vardecfn(16,5) = bigsigdk(3,3)/vardecfn(16,1);
      vardecfn(23,5) = bigsigdk(4,4)/vardecfn(23,1);
    end

    if k == 8
      vardecfn(3,5) = bigsigdk(1,1)/vardecfn(3,1);
      vardecfn(10,5) = bigsigdk(2,2)/vardecfn(10,1);
      vardecfn(17,5) = bigsigdk(3,3)/vardecfn(17,1);
      vardecfn(24,5) = bigsigdk(4,4)/vardecfn(24,1);
    end

    if k == 12
      vardecfn(4,5) = bigsigdk(1,1)/vardecfn(4,1);
      vardecfn(11,5) = bigsigdk(2,2)/vardecfn(11,1);
      vardecfn(18,5) = bigsigdk(3,3)/vardecfn(18,1);
      vardecfn(25,5) = bigsigdk(4,4)/vardecfn(25,1);
    end

    if k == 20
      vardecfn(5,5) = bigsigdk(1,1)/vardecfn(5,1);
      vardecfn(12,5) = bigsigdk(2,2)/vardecfn(12,1);
      vardecfn(19,5) = bigsigdk(3,3)/vardecfn(19,1);
      vardecfn(26,5) = bigsigdk(4,4)/vardecfn(26,1);
    end

    if k == 40
      vardecfn(6,5) = bigsigdk(1,1)/vardecfn(6,1);
      vardecfn(13,5) = bigsigdk(2,2)/vardecfn(13,1);
      vardecfn(20,5) = bigsigdk(3,3)/vardecfn(20,1);
      vardecfn(27,5) = bigsigdk(4,4)/vardecfn(27,1);
    end

    bigfx = bigfx*bigax;

  end    

  bigsigs1 = inv(eye(81)-kron(bigax,bigax))*bigbvbx(:);

  bigsigs = reshape(bigsigs1,9,9);

  bigsigd = bigcx*bigsigs*bigcx';

  vardecfn(7,5) = bigsigd(1,1)/vardecfn(7,1);
  vardecfn(14,5) = bigsigd(2,2)/vardecfn(14,1);
  vardecfn(21,5) = bigsigd(3,3)/vardecfn(21,1);
  vardecfn(28,5) = bigsigdk(4,4)/vardecfn(28,1);

% convert to percentages

  vardecfn = 100*vardecfn;