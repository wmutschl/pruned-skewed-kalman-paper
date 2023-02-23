% vardec.m: Uses the parameter estimates and standard errors found by
%             estse.m to compute variance decompositions for the New 
%             Keynesian model with markup and technology shocks.
%
%           Returns four 7x5 matrices: gvar, pvar, rvar, and xvar.  Each of
%             these matrices takes the form
%
%             _var = [ _vartot _vara _vare _varz _varr ],
%
%             where _vartot is a 7x1 vector of k-step ahead forecast
%             error variances in the variable _ (_=g for output, 
%             _=p for inflation, _=r for the nominal interest rate, or
%             _=x for the output gap), and where k = 1, 4, 8, 12, 20, 40,
%             and infinity. The variances are all expressed as percentages.
%             The 7x1 vectors _vara, _vare, _varz, and _varr are vectors of
%             the percentages of _vartot due to IS, markup, technology, and
%             policy shocks.
%
%           Also returns four 7x5 matrices: gvarse, pvarse, rvarse, and
%             xvarse containing standard errors for the corresponding elements
%             of gvar, pvar, rvar, zvar.
%
%           Formulas used to compute variance decompositions and standard
%             errors are contained in vardecfn.m.
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

% compute variance decompositions

  vdec = vardecfn(tstars);

  vdec1 = vdec(:);

% compute standard errors

  eee = 1e-6;

  epsmat = eee*eye(12);

  gradmat = zeros(140,12);

  for i = 1:12

    vdec2 = vardecfn(tstars+epsmat(:,i));

    vdec2 = vdec2(:);

    gradmat(:,i) = (vdec2-vdec1)/eee;

  end

  vdecse1 = gradmat*inv(hessmatf)*gradmat';

  vdecse2 = sqrt(diag(vdecse1));

  vdecse = reshape(vdecse2,28,5);

% create output matrices

  gvar = vdec(1:7,:);
  pvar = vdec(8:14,:);
  rvar = vdec(15:21,:);
  xvar = vdec(22:28,:);

  gvarse = vdecse(1:7,:);
  pvarse = vdecse(8:14,:);
  rvarse = vdecse(15:21,:);
  xvarse = vdecse(22:28,:);