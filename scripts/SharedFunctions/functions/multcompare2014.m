function [comparison,means,h,gnames] = multcompare(stats,varargin)
%MULTCOMPARE Perform a multiple comparison of means or other estimates
%   MULTCOMPARE performs a multiple comparison using one-way anova or
%   anocova results to determine which estimates (such as means,
%   slopes, or intercepts) are significantly different.
%
%   COMPARISON = MULTCOMPARE(STATS) performs a multiple comparison
%   using a STATS structure that is obtained as output from any of
%   the following functions:  anova1, anova2, anovan, aoctool,
%   kruskalwallis, friedman.  The return value COMPARISON is a matrix
%   with one row per comparison and six columns.  Columns 1-2 are the
%   indices of the two samples being compared.  Columns 3-5 are a lower
%   bound, estimate, and upper bound for their difference. Column 6 is the 
%   p-value for each individual comparison. 
%
%   COMPARISON = MULTCOMPARE(STATS, 'PARAM1',val1, 'PARAM2',val2,...)
%   specifies one or more of the following name/value pairs:
%   
%     'alpha'       Specifies the confidence level as 100*(1-ALPHA)%
%                   (default 0.05).
%     'display'     Either 'on' (the default) to display a graph of the
%                   estimates with comparison intervals around them, or
%                   'off' to omit the graph.
%     'ctype'       The type of critical value to use.  Choices are
%                   'tukey-kramer' (default), 'dunn-sidak', 'bonferroni',
%                   'scheffe'.  Enter two or more choices separated by
%                   spaces, to use the minimum of those critical values.
%     'dimension'   A vector specifying the dimension or dimensions over
%                   which the population marginal means are to be
%                   calculated.  Used only if STATS comes from anovan.
%                   The default is to compute over the first dimension
%                   associated with a categorical (non-continuous) factor.
%                   The value [1 3], for example, computes the population
%                   marginal mean for each combination of the first and
%                   third predictor values.
%     'estimate'    Estimate to compare.  Choices depend on the source of
%                   the stats structure:
%         anova1:  ignored, compare group means
%         anova2:  'column' (default) or 'row' means
%         anovan:  ignored, compare population marginal means
%         aoctool:  'slope', 'intercept', or 'pmm' (default is 'slope'
%                   for separate-slopes models, 'intercept' otherwise)
%         kruskalwallis:  ignored, compare average ranks of columns
%         friedman:  ignored, compare average ranks of columns
%
%   [COMPARISON,MEANS,H,GNAMES] = MULTCOMPARE(...) returns additional
%   outputs.  MEANS is a matrix with columns equal to the estimates
%   and their standard errors.  H is a handle to the figure containing
%   the graph.  GNAMES is a cell array with one row for each group,
%   containing the names of the groups.
%
%   The intervals shown in the graph are computed so that to a very close
%   approximation, two estimates being compared are significantly different
%   if their intervals are disjoint, and are not significantly different if
%   their intervals overlap.  (This is exact for multiple comparison
%   of means from anova1, if all means are based on the same sample size.)
%   You can click on any estimate to see which means are significantly
%   different from it.
%
%   Two additional CTYPE choices are available.  The 'hsd' option stands
%   for "honestly significant differences" and is the same as the
%   'tukey-kramer' option.  The 'lsd' option stands for "least significant
%   difference" and uses plain t-tests; it provides no protection against
%   the multiple comparison problem unless it follows a preliminary overall
%   test such as an F test.
%
%   MULTCOMPARE does not support multiple comparisons using anovan output
%   for a model that includes random or nested effects.  The calculations
%   for a random effects model produce a warning that all effects are
%   treated as fixed.  Nested models are not accepted.
%
%   Example:  Perform 1-way anova, and display group means with their names
%
%      load carsmall
%      [p,t,st] = anova1(MPG,Origin,'off');
%      [c,m,h,nms] = multcompare(st,'display','off');
%      [nms num2cell(m)]   
%
%   See also ANOVA1, ANOVA2, ANOVAN, AOCTOOL, FRIEDMAN, KRUSKALWALLIS.

%   Also supports older calling sequence:
%      [...] = MULTCOMPARE(STATS,ALPHA,DISPLAYOPT,CTYPE,ESTIMATE,DIM)
%
%   Reference: Y. Hochberg and A.C. Tamhane, "Multiple Comparison
%   Procedures," Wiley, New York, 1987.
%
%   The Tukey-Kramer critical value is the default.  This is known
%   to be the best choice for one-way anova comparisons, but the
%   conjecture that this is best for other comparisons is
%   unproven.  The Bonferroni and Scheffe critical values are
%   conservative in all cases.

%   Copyright 1993-2012 The MathWorks, Inc.
%   Copyright 1993-2013 The MathWorks, Inc.
%   Copyright 1993-2014 The MathWorks, Inc.


% Check inputs and assign default values
narginchk(1,Inf);

if isempty(varargin) || (~isempty(varargin{1}) && ischar(varargin{1}))
   okargs =   {'alpha' 'displayopt' 'ctype' 'estimate' 'dimension'};
   defaults = {0.05    'on'         ''      ''         []};
   [alpha,displayopt,ctype,estimate,dim] = ...
                internal.stats.parseArgs(okargs,defaults,varargin{:});
else
   % Old-style calling sequence with fixed argument positions
   if (nargin >= 2) && ~isempty(varargin{1})
      alpha = varargin{1};
   else
      alpha = 0.05;
   end
   if (nargin >= 3) && ~isempty(varargin{2})
      displayopt = varargin{2};
   else
      displayopt = 'on';
   end
   if (nargin>=4) && ~isempty(varargin{3})
      ctype = varargin{3};
   else
      ctype = '';
   end
   if (nargin>=5) && ~isempty(varargin{4})
      estimate = varargin{4};
   else
      estimate = '';
   end
   if (nargin>=6) && ~isempty(varargin{5})
      dim = varargin{5};
   else
      dim = [];
   end
end

if ~isstruct(stats) || ~isfield(stats,'source')
   error(message('stats:multcompare:BadStats'));
end
if (length(alpha)~=1 || ~isfinite(alpha) || alpha<=0 || alpha>=1)
   error(message('stats:multcompare:BadAlpha'));
end
if ~(isequal(displayopt, 'on') || isequal(displayopt, 'off'))
   error(message('stats:multcompare:BadDisplayOpt'));
end
dodisp = isequal(displayopt, 'on');

source = stats.source;

switch(source)
 case 'anova1'
   mname = 'means';
   gmeans = stats.means(:);
   gnames = stats.gnames;
   n = stats.n(:);
   df = stats.df;
   s = stats.s;
   ng = sum(n>0);
   if (df < 1)
      error(message('stats:multcompare:NotEnoughDataANOVA'));
   end
   
   gcov = diag((s^2)./n);
   [mn,se] = tValue(gmeans,gcov);
   t = mn./se;
   
   % Get Tukey-Kramer critical value
   if (isempty(ctype)), ctype = 'tukey-kramer'; end
   [crit,pval] = getcrit(ctype, alpha, df, ng, t);

 case 'anova2'
   docols = true;
   if (~isempty(estimate))
      estimate = internal.stats.getParamVal(estimate,{'row' 'column'},'ESTIMATE');
      docols = isequal(estimate,'column');
   end
   if (docols)
      gmeans = stats.colmeans(:);
      n = stats.coln(:);
      mname = 'column means';
   else
      gmeans = stats.rowmeans(:);
      n = stats.rown(:);
      mname = 'row means';
   end
   ng = length(gmeans);
   sigma = sqrt(stats.sigmasq);
   gnames = strjust(num2str((1:ng)'), 'left');
   df = stats.df;
   if (df < 1)
      error(message('stats:multcompare:NotEnoughDataANOVA'));
   end
   
   gcov = ((sigma^2)/n) * eye(ng);
   [mn,se] = tValue(gmeans,gcov);
   t = mn./se;
   
   % Get Tukey-Kramer critical value
   if (isempty(ctype)), ctype = 'tukey-kramer'; end
   [crit,pval] = getcrit(ctype, alpha, df, ng, t); 

   % This whole activity is a little strange if the model includes
   % interactions, especially if they are important.
   if (stats.inter && dodisp)     % model included an interaction term
      if (stats.pval < alpha)
         theString = getString(message('stats:multcompare:NoteSignifInteraction'));
      else
         theString = getString(message('stats:multcompare:NoteInsignifInteraction'));
      end
      theCell = textwrap({theString},80);
      fprintf('%s\n',theCell{:});
   end

 case 'anovan'
   mname = 'population marginal means';
   
   % We do not handle nested models
   try
       vnested = stats.vnested;
   catch me %#ok<NASGU>
       vnested = [];
   end
   if any(vnested(:))
       error(message('stats:multcompare:NoNesting'));
   end


   % Our calculations treat all effects as fixed
   if ~isempty(stats.ems)
      warning(message('stats:multcompare:IgnoringRandomEffects'))
   end
   
   nvars = length(stats.nlevels);
   P0 = stats.nullproject;
   
   % Make sure DIM is a scalar or vector of factor numbers.
   if isempty(dim)
       dim = find(stats.nlevels>1,1);
   end
   dim = dim(:);
   if isempty(dim) || any(dim<1 | dim>nvars | dim~=round(dim))
      error(message('stats:multcompare:BadDim', nvars));
   end
   dim = sort(dim);
   dim(diff(dim)==0) = [];
   if any(stats.nlevels(dim)<2)
      error(message('stats:multcompare:DimSpecifiesZeroDFFactor'));
   end
   
   % Create all combinations of the specified factors
   try
       continuous = stats.continuous;
   catch me %#ok<NASGU>
       continuous = zeros(nvars,1);
   end
   ffdesign = fullfact(stats.nlevels(dim));
   nrows = size(ffdesign,1);
   
   % Create a design matrix for these combinations
   dums = cell(nvars, 1);
   for j=1:length(dim)
      dj = dim(j);
      dums{dj} = idummy(ffdesign(:,j),3);  % dummy variables for each factor
   end
   
   % Create a vector of average values for remaining factors
   for j=1:nvars
      if isempty(dums{j});
          if continuous(j)
             dums{j} = stats.vmeans(j) * ones(nrows,1);
          else
             nlev = stats.nlevels(j);
             dums{j} = (1/nlev) * ones(nrows,nlev);
          end
      end
   end

   % Fill in x columns for each term
   termcols = stats.termcols(:);
   termstart = cumsum([1; termcols]);
   terms = [zeros(1,nvars); stats.terms];
   ncols = sum(termcols);
   x = zeros(size(ffdesign,1), ncols);
   x(:,1) = 1;
   for j=1:length(termcols)
      tm = terms(j,:);
      t0 = termstart(j);
      t1 = termstart(j) + termcols(j) - 1;
      if all(tm==0)
         x(:,t0:t1) = 1;
      else
         x0 = [];
         for k=nvars:-1:1
            if tm(k)
               x0 = termcross(x0,dums{k});
            end
         end
         x(:,t0:t1) = x0;
      end
   end

   % Compute estimates and their standard errors
   gmeans = x * stats.coeffs;
   xproj = (x*P0)';
   tmp = stats.Rtr \ xproj;
   if (stats.dfe == 0)
      mse = NaN;
   else
      mse = max(stats.mse,0);
   end
   gcov = mse * (tmp' * tmp);
   
   % Find non-estimable means and set them to NaN
   Xrows = stats.rowbasis';           % row basis of original X matrix
   bb = Xrows \ (x');                 % fit rows of x to row basis
   xres = Xrows * bb - x';            % get residuals
   xres = sum(abs(xres));             % sum of absolute residuals
   cutoff = sqrt(eps(class(xres))) * size(xres,2); % cutoff for large residuals
   gmeans(xres > cutoff) = NaN;       % not in row space of original X
   
   [mn,se] = tValue(gmeans,gcov);
   t = mn./se;
   
   % Get Tukey-Kramer critical value
   if (isempty(ctype)), ctype = 'tukey-kramer'; end
   [crit,pval] = getcrit(ctype, alpha, stats.dfe, length(gmeans),t);

   % Get names for each group
   ngroups = size(ffdesign,1);
   gnames = cell(ngroups,1);
   allnames = stats.grpnames;
   varnames = stats.varnames;
   for j=1:ngroups
      v1 = dim(1);
      vals1 = allnames{v1};
      nm = sprintf('%s=%s',varnames{v1},vals1{ffdesign(j,1)});
      for i=2:size(ffdesign,2)
         v2 = dim(i);
         vals2 = allnames{v2};
         nm = sprintf('%s,%s=%s',nm,varnames{v2},vals2{ffdesign(j,i)});
      end
      gnames{j} = nm;
   end

 case 'aoctool'
   model = stats.model;
   if (model==1 || model==3)
      error(message('stats:multcompare:NoMultipleParameters'));
   end
   gnames = stats.gnames;
   n = stats.n(:);
   ng = length(n);
   df = stats.df;
   if (df < 1)
      error(message('stats:multcompare:NotEnoughDataAOC'));
   end

   % Get either slope or intercept estimates and covariances
   if (isempty(estimate))
      if (model == 5)
         estimate = 'slope';
      else
         estimate = 'intercept';
      end
   else
      estimate = internal.stats.getParamVal(estimate,{'slope' 'intercept' 'pmm'},'ESTIMATE');
   end
   switch(estimate)
    case 'slope'
      if (~isfield(stats, 'slopes'))
         error(message('stats:multcompare:BadStatsNoSlope'));
      end
      gmeans = stats.slopes;
      gcov = stats.slopecov;
      mname = 'slopes';
    case 'intercept'
      if (~isfield(stats, 'intercepts'))
         error(message('stats:multcompare:BadStatsNoIntercept'));
      end
      gmeans = stats.intercepts;
      gcov = stats.intercov;
      mname = 'intercepts';
    case 'pmm'
      gmeans = stats.pmm;
      gcov = stats.pmmcov;
      mname = 'population marginal means';
   end

   if (any(any(isinf(gcov))))
      error(message('stats:multcompare:InfiniteVariance', mname));
   end
   
   [mn,se] = tValue(gmeans,gcov);
   t = mn./se;
   
   % Get Tukey-Kramer critical value
   if (isempty(ctype)), ctype = 'tukey-kramer'; end
   [crit,pval] = getcrit(ctype, alpha, df, ng, t);

 case 'kruskalwallis'
   gmeans = stats.meanranks(:);
   gnames = stats.gnames;
   n = stats.n(:);
   sumt = stats.sumt;
   ng = length(n);
   N = sum(n);
   mname = 'mean ranks';

   gcov = diag(((N*(N+1)/12) - (sumt/(12*(N-1)))) ./ n);
   [mn,se] = tValue(gmeans,gcov);
   t = mn./se;
   
   % Get critical value; H&T recommend the Tukey-Kramer value
   if (isempty(ctype)), ctype = 'tukey-kramer'; end
   [crit,pval] = getcrit(ctype, alpha, Inf, ng, t);
   
   % Note that the intervals in M can be used for testing but not
   % for simultaneous confidence intervals.  See H&T, p. 249.
   if (dodisp)
      disp(getString(message('stats:multcompare:NoteNotSimul')));
   end

 case 'friedman'
   gmeans = stats.meanranks(:);
   n = stats.n;
   ng = length(gmeans);
   sigma = stats.sigma;
   mname = 'mean column ranks';
   gnames = strjust(num2str((1:ng)'), 'left');

   gcov = ((sigma^2) / n) * eye(ng);
   [mn,se] = tValue(gmeans,gcov);
   t = mn./se;
   
   % Get critical value; H&T recommend the Tukey-Kramer value
   if (isempty(ctype)), ctype = 'tukey-kramer'; end
   [crit,pval] = getcrit(ctype, alpha, Inf, ng, t);

   % Note that the intervals in M can be used for testing but not
   % for simultaneous confidence intervals.  See H&T, p. 249.
   if (dodisp)
      disp(getString(message('stats:multcompare:NoteNotSimul')));
   end

 otherwise
   error(message('stats:multcompare:BadStats'));
end


% Create output matrix showing tests for all pairwise comparisons
% and graph that approximates these tests.
[M,MM,hh] = makeM(gmeans, gcov, mn, se, crit, gnames, mname, dodisp, pval);

comparison = M;
if (nargout>1), means = MM; end
if (nargout>2), h = hh; end

% -----------------------------------------------
function [crit,pval] = getcrit(ctype, alpha, df, ng, t)
% Get the minimum of the specified critical values
crit = Inf;
[onetype,ctype] = strtok(ctype);

while(~isempty(onetype))
   if (length(onetype) == 1)
      switch onetype
       case 't', onetype = 'tukey-kramer';
       case 'd', onetype = 'dunn-sidak';
       case 'b', onetype = 'bonferroni';
       case 's', onetype = 'scheffe';
       case 'h', onetype = 'tukey-kramer';
       case 'l', onetype = 'lsd';
      end
   end
   if (isequal(onetype, 'hsd')), onetype = 'tukey-kramer'; end
   
   switch onetype
    case 'tukey-kramer' % or hsd
     crit1 = stdrinv(1-alpha, df, ng) / sqrt(2);
     
     % The T-K algorithm is inaccurate for small alpha, so compute
     % an upper bound for it and make sure it's in range.
     ub = getcrit('dunn-sidak', alpha, df, ng, t);
     if (crit1 > ub), crit1 = ub; end
     pval = 0*t;
     for j=1:numel(t)
         pval(j) = 1 - stdrcdf(sqrt(2)*abs(t(j)),df,ng);
     end

    case 'dunn-sidak'
     kstar = nchoosek(ng, 2);
     alf = 1-(1-alpha).^(1/kstar);
     if (isinf(df))
        crit1 = norminv(1-alf/2);
     else
        crit1 = tinv(1-alf/2, df);
     end
     pval = 1 - (1-2*tcdf(-abs(t),df)).^kstar;

    case 'bonferroni'
     kstar = nchoosek(ng, 2);
     if (isinf(df))
        crit1 = norminv(1 - alpha / (2*kstar));
     else
        crit1 = tinv(1 - alpha / (2*kstar), df);
     end
     pval = 2*kstar*tcdf(-abs(t),df);

    case 'lsd'
     if (isinf(df))
        crit1 = norminv(1 - alpha / 2);
     else
        crit1 = tinv(1 - alpha / 2, df);
     end
     pval = 2*tcdf(-abs(t),df);

    case 'scheffe'
     if (isinf(df))
        tmp = chi2inv(1-alpha, ng-1) / (ng-1);
     else
        tmp = finv(1-alpha, ng-1, df);
     end
     crit1 = sqrt((ng-1) * tmp);
     pval = fcdf((t.^2)/(ng-1),ng-1,df,'upper');
     
    otherwise
     error(message('stats:multcompare:BadCType', ctype));
   end

   pval(pval>1) = 1;
   if (~isnan(crit1)), crit = min(crit, crit1); end
   [onetype,ctype] = strtok(ctype);
end

% -----------------------------------------------
function [M,MM,hh] = makeM(gmeans, gcov, mn, se, crit, gnames, mname, dodisp, pval)
% Create matrix to test differences, matrix of means, graph to display test

ng = length(gmeans);
MM = zeros(ng,2);
MM(:,1) = gmeans(:);
MM(:,2) = sqrt(diag(gcov));
MM(isnan(MM(:,1)),2) = NaN;

M = nchoosek(1:ng, 2);      % all pairs of group numbers
M(1,5) = 0;                 % expand M to proper size
g1 = M(:,1);
g2 = M(:,2);
M(:,4) = mn;
i12 = sub2ind(size(gcov), g1, g2);
delta = crit * se;
M(:,3) = M(:,4) - delta;
M(:,5) = M(:,4) + delta;
M(:,6) = pval;

% If requested, make a graph that approximates these tests
if (dodisp)
   % Find W values according to H&T (3.32, p. 98)
   d = zeros(ng, ng);
   d(i12) = se;
   sum1 = sum(sum(d));
   d = d + d';
   sum2 = sum(d);
   if (ng > 2)
      w = ((ng-1) * sum2 - sum1) ./ ((ng-1)*(ng-2));
   else
      w = sum1 * ones(2, 1) / 2;
   end
   halfwidth = crit * w(:);
   hh = meansgraph(gmeans, gmeans-halfwidth, gmeans+halfwidth, ...
                   gnames, mname);
   set(hh, 'Name', getString(message('stats:multcompare:MultcompareFigureTitleString', mname)));
else
   hh = [];
end

% -----------------------------------------------
function [mn,se]= tValue(gmeans,gcov)
% Make sure NaN groups don't affect other results
t = isnan(gmeans);
if any(t)
    gcov(t,:) = 0;
    gcov(:,t) = 0;
end
ng = length(gmeans);
M = nchoosek(1:ng, 2);      % all pairs of group numbers
g1 = M(:,1);
g2 = M(:,2);
mn = gmeans(g1) - gmeans(g2);
i12 = sub2ind(size(gcov), g1, g2);
gvar = diag(gcov);
se = sqrt(gvar(g1) + gvar(g2) - 2 * gcov(i12));              


function xout = stdrcdf(q, v, r, upper)
%STDRCDF Compute c.d.f. for Studentized Range statistic
%   F = STDRCDF(Q,V,R) is the cumulative distribution function for the
%   Studentized range statistic for R samples and V degrees of
%   freedom, evaluated at Q.
%
%   G = STDRCDF(Q,V,R,'upper') is the upper tail probability,
%   G=1-F.  This version computes the upper tail probability
%   directly (not by subtracting it from 1), and is likely to be
%   more accurate if Q is large and therefore F is close to 1.

%   Copyright 1993-2007 The MathWorks, Inc. 
%   $Revision: 1.1.8.2 $  $Date: 2010/10/08 17:29:44 $

% Based on Fortran program from statlib, http://lib.stat.cmu.edu
% Algorithm AS 190  Appl. Statist. (1983) Vol.32, No. 2
% Incorporates corrections from Appl. Statist. (1985) Vol.34 (1)
% Vectorized and simplified for MATLAB.  Added 'upper' option.

if (length(q)>1 | length(v)>1 | length(r)>1),
   error(message('stats:stdrcdf:NotScalar')); % for now
end
[err,q,v,r] = distchck(3,q,v,r);
if (err > 0)
   error(message('stats:stdrcdf:InputSizeMismatch'));
end
uppertail = 0;
if (nargin>3)
   if ~(  isequal(upper,'u') | isequal(upper,'upper') ...
        | isequal(upper,'l') | isequal(upper,'lower'))
      error(message('stats:stdrcdf:BadUpper'));
   end
   uppertail = isequal(upper,'u') | isequal(upper,'upper');
end

% Accuracy can be increased by use of a finer grid.  Increase
% jmax, kmax and 1/step proportionally.
jmax = 15;          % controls maximum number of steps
kmax = 15;          % controls maximum number of steps
step = 0.45;        % node spacing
vmax = 120;         % max d.f. for integration over chi-square

% Handle illegal or trivial values first.
xout = zeros(size(q));
if (length(xout) == 0), return; end   
ok = (v>0) & (v==round(v)) & (r>1) & (r==round(r));
xout(~ok) = NaN;
ok = ok & (q > 0);
v = v(ok);
q = q(ok);
r = r(ok);
if (length(v) == 0), return; end
xx = zeros(size(v));

% Compute constants, locate midpoint, adjust steps.
g = step ./ (r .^ 0.2);
if (v > vmax)
   c = log(r .* g ./ sqrt(2*pi));
else
   h = step ./ sqrt(v);
   v2 = v * 0.5;
   c = log(sqrt(2/pi) .* r .* g .* h) - v2 + v2.*log(v2) - gammaln(v2);

   j=(-jmax:jmax)';
   hj = h * j;
   ehj = exp(hj);
   qw = q .* ehj;
   vw = v .* (hj + 0.5 * (1 - ehj .^2));
   C = ones(1,2*kmax+1);         % index to duplicate columns
   R = ones(1,2*jmax+1);         % index to duplicate rows
end

% Compute integral by summing the integrand over a
% two-dimensional grid centered approximately near its maximum.
gk = (0.5 * log(r)) + g * (-kmax:kmax);
w0 = c - 0.5 * gk .^ 2;
pz = normcdf(-gk);
if (~uppertail)
   % For regular cdf, use integrand as in AS 190.
   if (v > vmax)
      % don't integrate over chi-square
      x = normcdf(q - gk) - pz;
      xx = sum(exp(w0) .* (x .^ (r-1)));
   else
      % integrate over chi-square
      x = normcdf(qw(:,C) - gk(R,:)) - pz(R,:);
      xx = sum(sum(exp(w0(R,:) + vw(:,C)) .* (x .^ (r-1))));
   end
else
   % To compute the upper tail probability, we need an integrand that
   % contains the normal probability of a region consisting of a
   % hyper-quadrant minus a rectangular region at the origin of the
   % hyperquadrant.
   if (v > vmax)           % for large d.f., don't integrate over chi-square
      xhq   = (1 - pz) .^ (r-1);
      xrect = (normcdf(q - gk) - pz) .^ (r-1);
      xx = sum(exp(w0) .* (xhq - xrect));
   else                    % for typical cases, integrate over chi-square
      xhq   = (1 - pz) .^ (r-1);
      xrect = (normcdf(qw(:,C) - gk(R,:)) - pz(R,:)) .^ (r-1);
      xx = sum(sum(exp(w0(R,:) + vw(:,C)) .* (xhq(R,:) - xrect)));
   end
end

xout(ok) = xx;



function d = idummy(x, method)
%DUMMY  Creates a matrix of dummy variables for a discrete variable
%   D=IDUMMY(X,METHOD) creates an array D of dummy variables for the
%   grouping variable I (integers 1,...,g), using the method specified:
%
%   method = 1:   0/-1/1 coding, full rank
%   method = 2:   0/1 coding, full rank
%   method = 3:   0/1 coding, overdetermined

%   Copyright 1993-2005 The MathWorks, Inc. 
%   $Revision: 1.1.8.1 $  $Date: 2010/03/16 00:29:42 $

if (nargin < 2)
   method = 1;
end

n = length(x);
g = max(x);
ncols = g - (method ~= 3);
d = zeros(n, ncols);

if (g > 1 || method==3)
   % Fill in -1 for the first level
   if (method == 1)
      d((x == 1),:) = -1;
   end
   
   % Fill in 1 in the appropriate column for other levels
   m3 = (method == 3);
   for j=(2-m3):g
      d((x == j),j-1+m3) = 1;
   end
end



function ab = termcross(a,b)
%TERMCROSS Multiply dummy variables for two terms to get interaction

%   Copyright 1993-2004 The MathWorks, Inc. 
%   $Revision: 1.1.8.1 $  $Date: 2010/03/16 00:30:36 $
if (isempty(a)), ab = b; return, end
if (isempty(b)), ab = a; return, end

na = size(a,2);
nb = size(b,2);
acols = repmat((1:na), 1, nb);
bcols = reshape(repmat((1:nb), na, 1), 1, na*nb);
ab = a(:,acols) .* b(:,bcols);




function x = stdrinv(p, v, r)
%STDRINV Compute inverse c.d.f. for Studentized Range statistic
%   STDRINV(P,V,R) is the inverse cumulative distribution function for
%   the Studentized range statistic for R samples and V degrees of
%   freedom, evaluated at P.

%   Copyright 1993-2013 The MathWorks, Inc. 


% Based on Fortran program from statlib, http://lib.stat.cmu.edu
% Algorithm AS 190  Appl. Statist. (1983) Vol.32, No. 2
% Incorporates corrections from Appl. Statist. (1985) Vol.34 (1)

if (length(p)>1 || length(v)>1 || length(r)>1),
   error(message('stats:stdrinv:NotScalar')); % for now
end

[err,p,v,r] = distchck(3,p,v,r);
if (err > 0)
   error(message('stats:stdrinv:InputSizeMismatch'));
end

% Handle illegal or trivial values first.
x = zeros(size(p));
if isempty(x), return; end
ok = (v>0) & (v==round(v)) & (r>1) & (r==round(r) & (p<1));
x(~ok) = NaN;
ok = ok & (p>0);
v = v(ok);
p = p(ok);
r = r(ok);
if isempty(v), return; end

% Define constants
jmax = 20;
pcut = 0.00001;
tiny = 0.000001;
upper = (p > .99);
if (upper)
   uppertail = 'u';
   p0 = 1-p;
else
   uppertail = 'l';
   p0 = p;
end

% Obtain initial values
q1 = qtrng0(p, v, r);
% p1 = internal.stats.stdrcdf(q1, v, r, uppertail);
p1 = stdrcdf(q1, v, r, uppertail);
xx = q1;
if (abs(p1-p0) >= pcut*p0)
   if (p1 > p0), p2 = max(.75*p0, p0-.75*(p1-p0)); end
   if (p1 < p0), p2 = p0 + (p0 - p1) .* (1 - p0) ./ (1 - p1) * 0.75; end
   if (upper)
      q2 = qtrng0(1-p2, v, r);
   else
      q2 = qtrng0(p2, v, r);
   end

   % Refine approximation
   for j=2:jmax
%       p2 = internal.stats.stdrcdf(q2, v, r, uppertail);
      p2 = stdrcdf(q2, v, r, uppertail);
      e1 = p1 - p0;
      e2 = p2 - p0;
      d = e2 - e1;
      xx = (q1 + q2) / 2;
      if (abs(d) > tiny*p0)
         xx = (e2 .* q1 - e1 .* q2) ./ d;
      end
      if (abs(e1) >= abs(e2))
         q1 = q2;
         p1 = p2;
      end
      if (abs(p1 - p0) < pcut*p0), break; end
	   q2 = xx;
   end
end
   
x(ok) = xx;

% ---------------------------------
function x = qtrng0(p, v, r)
% Algorithm AS 190.2  Appl. Statist. (1983) Vol.32, No.2
% Calculates an initial quantile p for a studentized range
% distribution having v degrees of freedom and r samples
% for probability p, p.gt.0.80 .and. p.lt.0.995.

t=norminv(0.5 + 0.5 .* p);
if (v < 120), t = t + 0.25 * (t.^3 + t) ./ v; end
q = 0.8843 - 0.2368 .* t;
if (v < 120), q = q - (1.214./v) + (1.208.*t./v); end
x = t .* (q .* log(r-1) + 1.4142);



% function x = stdrinv(p, v, r)
% %STDRINV
% % calls matlab's stdrinv function
% %
% % Example:
% %   if you don't have matlab statitics toolbox then edit this
% %   to change to your own strcdf function
% % $Id: stdrinv.m,v 1.3 2006/12/26 22:53:21 Mike Exp $
% % Copyright 2006 Mike Boedigheimer
% % Amgen Inc.
% % Department of Computational Biology
% % mboedigh@amgen.com
% % 
% x=dfswitchyard('stdrinv',p,v,r);
% 
% function xout = stdrcdf(q, v, r, upper)
% % STRCDF calls matlab's stdrcdf function
% %
% % Example:
% %   if you don't have matlab statitics toolbox then edit this
% %   to change to your own strcdf function
% % $Id: stdrcdf.m,v 1.3 2006/12/26 22:53:21 Mike Exp $
% % Copyright 2006 Mike Boedigheimer
% % Amgen Inc.
% % Department of Computational Biology
% % mboedigh@amgen.com
% % 
% upper = 'lower';
% xout=dfswitchyard('stdrcdf',q,v,r,upper);

function p = fcdf(x,v1,v2,uflag)
%FCDF   F cumulative distribution function.
%   P = FCDF(X,V1,V2) returns the F cumulative distribution function
%   with V1 and V2 degrees of freedom at the values in X.
%
%   The size of P is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.    
%
%   P = FCDF(X,V1,V2,'upper') returns the upper tail probability of the
%   F distribution with V1 and V2 degrees of freedom at the values in X.
%
%   See also FINV, FPDF, FRND, FSTAT, CDF.

%   References:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.6.

%   Copyright 1993-2013 The MathWorks, Inc. 


if nargin < 3, 
    error(message('stats:fcdf:TooFewInputs')); 
end

[errorcode, x, v1, v2] = distchck(3,x,v1,v2);

if errorcode > 0
    error(message('stats:fcdf:InputSizeMismatch'));
end

if nargin>=4
    if ~strcmpi(uflag,'upper')
        error(message('stats:cdf:UpperTailProblem'));
    end
    k = ~isnan(x);
    x(k) = 1./max(0,x(k));
    a=v1;
    b=v2;
    v1=b;
    v2=a;
end

% Initialize P to zero.
if isa(x,'single') || isa(v1,'single') || isa(v2,'single')
   p = zeros(size(x),'single');
else
   p = zeros(size(x));
end

t = (v1 <= 0 | v2 <= 0 | isnan(x) | isnan(v1) | isnan(v2));
p(t) = NaN;
s = (x==Inf) & ~t;
if any(s(:))
   p(s) = 1;
   t = t | s;
end

% Compute P when X > 0.
k = find(x > 0 & ~t & isfinite(v1) & isfinite(v2));
if any(k)
    k1 = (v2(k) <= x(k).*v1(k));
    % use A&S formula 26.6.2 to relate to incomplete beta function 
    % Also use 26.5.2 to avoid cancellation by subtracting from 1
    if any(k1)
        kk = k(k1);
        xx = v2(kk)./(v2(kk)+x(kk).*v1(kk));
        p(kk) = betainc(xx, v2(kk)/2, v1(kk)/2,'upper');
    end
    if any(~k1)
        kk = k(~k1);
        num = v1(kk).*x(kk);
        xx = num ./ (num+v2(kk));
        p(kk) = betainc(xx, v1(kk)/2, v2(kk)/2,'lower');
    end
end

if any(~isfinite(v1(:)) | ~isfinite(v2(:)))
   k = find(x > 0 & ~t & isfinite(v1) & ~isfinite(v2) & v2>0);
   if any(k)       
      p(k) = gammainc(v1(k).*x(k)./2, v1(k)./2, 'lower'); % chi2cdf(v1(k).*x(k),v1(k))
   end
   k = find(x > 0 & ~t & ~isfinite(v1) & v1>0 & isfinite(v2));
   if any(k)
      p(k) = gammainc(v2(k)./x(k)./2, v2(k)./2, 'upper'); % 1 - chi2cdf(v2(k)./x(k),v2(k))
   end
   k = find(x > 0 & ~t & ~isfinite(v1) & v1>0 & ~isfinite(v2) & v2>0);
   if any(k)
       if nargin>=4 && x(k)==1
           p(k) = 0;
       else
           p(k) = (x(k)>=1);
       end
   end
end  

