function [ p, table, stats, terms ] = myanovan ( data, group, varargin )
%ANOVAN N-way analysis of variance (ANOVA).
%   P=ANOVAN(Y,GROUP) performs multi-way (n-way) anova on the vector Y
%   grouped by entries in the cell array GROUP.  Each cell of GROUP must
%   contain a grouping variable that can be a categorical variable, numeric
%   vector, character matrix, or single-column cell array of strings.
%   Each grouping variable must have the same number of items as Y.  The
%   fitted anova model includes the main effects of each grouping variable.
%   All grouping variables are treated as fixed effects by default.  The
%   result P is a vector of p-values, one per term.
%
%   P=ANOVAN(Y,GROUP,'PARAM1',val1,'PARAM2',val2,...) specifies one or
%   more of the following name/value pairs:
%
%     Parameter    Value
%     'alpha'      A value between 0 and 1 requesting 100*(1-alpha)%
%                  confidence bounds (default 0.05 for 95% confidence)
%     'continuous' A vector of indices indicating which grouping variables
%                  should be treated as continuous predictors rather than
%                  as categorical predictors
%     'display'    Either 'on' (the default) to display an anova table
%                  or 'off' to omit the display
%     'nested'     A matrix M of 0's and 1's specifying the nesting
%                  relationships among the grouping variables.  M(i,j)
%                  is 1 if variable i is nested in variable j.
%     'random'     A vector of indices indicating which grouping variables
%                  are random effects (all are fixed by default)
%     'sstype'     The type of sum of squares 1, 2, 3, or 'h' (default=3)
%     'varnames'   Grouping variables names in a character matrix or
%                  a cell array of strings, one per grouping variable
%                  (default names are 'X1', 'X2', ...)
%
%     'model'      The model to use, specified as one of the following:
%
%        'linear' to use only main effects of all factors (default)
%        'interaction' for main effects plus two-factor interactions
%        'full' to include interactions of all levels
%        an integer representing the maximum interaction order, for example
%           3 means main effects plus two- and three-factor interactions
%        a matrix of term definitions as accepted by the X2FX function,
%           but all entries must be 0 or 1 (no higher powers)
%
%   [P,T,STATS,TERMS]=ANOVAN(...) also returns a cell array T containing
%   the anova table, a structure STATS containing a variety of statistics,
%   and a matrix TERMS describing the terms used (suitable for use as the
%   MODEL input argument if you run ANOVAN again).
%
%   For models without random effects, the anova table T contains columns
%   for terms, sum of squares, degrees of freedom, an indication of whether
%   the term is singular, mean square, F statistic, and P value.  For models
%   with random effects there are additional columns showing the term type
%   (fixed or random), expected mean square, denominator mean square for F,
%   denominator degrees of freedom for F, denominator definition, variance
%   component estimate, lower bound for variance, and upper bound for
%   variance.
%
%   For the 'sstype' parameter, values 1-3 produce the usual type 1, type 2,
%   or type 3 sums of squares.  The value 'h' produces sums of squares for
%   a hierarchical model that is similar to type 2, but with continuous as
%   well as categorical factors used to determine the hierarchy of terms.
%
%   The STATS structure contains the fields listed below, in addition
%   to a number of other fields required for doing multiple comparisons
%   using the MULTCOMPARE function:
%      coeffs      estimated coefficients
%      coeffnames  name of term for each coefficient
%      vars        matrix of grouping variable values for each term
%      resid       residuals from the fitted model
%
%   These fields exist if there are random effects:
%      ems         expected mean squares
%      denom       denominator definition
%      rtnames     names of random terms
%      varest      variance component estimates (one per random term)
%      varci       confidence intervals for variance components
%
%   Example:
%    load carbig
%    [p, atab] = anovan(MPG, {Cylinders Origin Model_Year}, ...
%                       'model',2, 'sstype',2, ...
%                       'varnames',strvcat('Cyl', 'Org', 'Yr'))
%    This example performs three-way anova on MPG using the factors
%    Cylinders, Origin, and Year.  The model will have all two-factor
%    interactions but not the three-factor interaction.  Sum of
%    squares are Type 2.
%
%   See also MULTCOMPARE, ANOVA1, ANOVA2, MANOVA1.

%   Also supports older calling sequence:
%      P=ANOVAN(Y,GROUP,MODEL,SSTYPE,VARNAME,DISPLAYOPT)
%   In addition to the matrix form of the model specification, the function
%   allows MODEL to be a vector V of integers.  In this compressed form
%   each element describes a term, so the Ith term includes the Jth grouping
%   variable if BITGET(V(I),J)==1
%
%
%   References:
%      Dunn, O.J., and V.A. Clark (1974), "Applied Statistics:
%         Analysis of Variance and Regression," Wiley, New York.
%      Goodnight, J.H., and F.M. Speed (1978), "Computing Expected Mean
%         Squares," SAS Institute, Cary, NC.
%      Milliken, G.A., and D.E. Johnson, "Analysis of Messy Data,"
%         Chapman & Hall, New York.
%      Seber, G.A.F. (1977), "Linear Regression Analysis," Wiley,
%         New York.

%   Copyright 1993-2012 The MathWorks, Inc.
%   $Revision: 1.1.8.6 $  $Date: 2012/04/02 22:16:41 $

% Checks that at least there are to input arguments.
narginchk ( 2, Inf );

% If data is a vector transforms it to a column vector.
if numel ( data ) == length ( data )
    data = data (:);
end

% Gets the number of cases and the shape of the comparison matrix.
mshape = size ( data );
mshape (1) = 1;

% Reshapes the comparison matrix to 2D.
data   = data ( :, : );
ncomps = size ( data, 2 );

if ~iscell ( group ) && ~isa ( group, 'categorical' )
    if isnumeric ( group ) && size ( group, 1 ) == ncases
        group = num2cell ( group, 1 );
    else
        error ( message ( 'stats:anovan:GroupNotCellOrGroupingMatrix' ) )
    end
end
group  = group (:);
ngroup = length ( group );

% Initializes the input arguments.
model      = 1;
sstype     = 3;
varnames   = '';
displayopt = 'on';
alpha      = 0.05;
randomvar  = false ( ngroup, 1 );
continuous = [];
vnested    = [];
dohier     = false;

% Enables the legacy style, if required.
if nargin > 2 && nargin <= 6 && ( ~ischar ( varargin {1} ) || isempty ( varargin {1} ) || ismember ( varargin {1}, { 'l' 'i' 'f' 'linear' 'interaction' 'full' } ) )
    legacy = true;
else
    legacy = false;
end

% If legacy gets the input variables in in fixed positions.
if legacy
    if nargin >= 3, model      = varargin {1}; end
    if nargin >= 4, sstype     = varargin {2}; end
    if nargin >= 5, varnames   = varargin {3}; end
    if nargin >= 6, displayopt = varargin {4}; end
    
% Otherwise gets the variables as key-value pairs.
else
    
    % Gets the arguments provided.
    par1 = 1 + find ( strcmp ( varargin, 'model' ) );
    par2 = 1 + find ( strcmp ( varargin, 'sstype' ) );
    par3 = 1 + find ( strcmp ( varargin, 'varnames' ) );
    par4 = 1 + find ( strcmp ( varargin, 'display' ) );
    par5 = 1 + find ( strcmp ( varargin, 'alpha' ) );
    par6 = 1 + find ( strcmp ( varargin, 'random' ) );
    par7 = 1 + find ( strcmp ( varargin, 'continuous' ) );
    par8 = 1 + find ( strcmp ( varargin, 'nested' ) );
    
    % Updates the arguments.
    if ~isempty ( par1 ), model      = varargin { par1 }; end
    if ~isempty ( par2 ), sstype     = varargin { par2 }; end
    if ~isempty ( par3 ), varnames   = varargin { par3 }; end
    if ~isempty ( par4 ), displayopt = varargin { par4 }; end
    if ~isempty ( par5 ), alpha      = varargin { par5 }; end
    if ~isempty ( par6 ), randomvar  = varargin { par6 }; end
    if ~isempty ( par7 ), continuous = varargin { par7 }; end
    if ~isempty ( par8 ), vnested    = varargin { par8 }; end
end

% Check optional arguments
[varnames,termlist,randomvar,sstype,continuous,vnested] = argcheck(ngroup,...
    varnames,model,randomvar,sstype,alpha,continuous,vnested);

% Sets the hierarchical sum of squares to type 2.
if isequal ( lower ( sstype ), 'h' )
    sstype = 2;
    dohier = true;
end

% If there is any random variable computes expected mean of squares.
doems = any ( randomvar );

% Rewrites the nested variable names.
if any ( vnested (:) )
    varnames = makenestednames ( varnames, vnested );
end


% STEP 1:  Remove NaN's and prepare grouping variables.
% Also, make sure all groups are rows, and create group index and name arrays.
[ data, allgrps, nanrow, varinfo ] = removenans ( data, group, continuous );
ncases = size ( data, 1 );


% STEP 2:  Create dummy variable arrays for each grouping variable.
% STEP 2a: For type 3 ss, create constraints for each grouping variable.
varinfo.varnames = varnames;
varinfo = makedummyvars ( continuous, allgrps, vnested, varinfo );


% STEP 3:  Create dummy variable arrays for each term in the model.
nterms = size ( termlist, 1 );
fulltermlist = showtermnesting ( termlist, vnested );
tnested = gettermnesting ( fulltermlist, vnested, continuous );

[ ~, sindex ] = sortrows ( fulltermlist );
terminfo = makedummyterms ( sindex, termlist, varinfo, continuous, tnested );


% STEP 4:  Create the full design matrix
[ dmat, cmat, termname ] = makedesignmat ( terminfo, ncases );


% STEP 5:  Fit the full model and compute the residual sum of squares
% Remvoes the mean in two passes (to improve accuracy).
mu   =  mean ( data, 1 );
data = bsxfun ( @minus, data, mu );
mu2  = mean ( data, 1 );
mu   = mu + mu2;
data = bsxfun ( @minus, data, mu2 );
sst  = sum ( data .^ 2, 1 );

if nargout >= 3
    % Do initial fit requesting stats structure, then convert from full-rank
    % to overdetermined form, then find a label for each coefficient
    [ ssx, dfx, dmat2, data, stats ] = dofit ( dmat, data, cmat, sst, mu );
    [ codes, names ] = convertstats ( stats, terminfo, varinfo, continuous );
else
    [ ssx, dfx, dmat2, data ] = dofit ( dmat, data, cmat, sst, mu );
end

sse  = sst - ssx;
dfe  = ncases - dfx;


% STEP 6:  Determine which models to compare for testing each term
ssw  = -ones ( nterms, ncomps );      % sum of squares with this term
sswo = ssw;                   % sum of squares without this term
dfw  = ssw;                   % residual d.f. with this term
dfwo = ssw;                   % residual d.f. without this term

switch sstype
    case 1
        modw = tril ( true ( nterms ) ); % get model with this term
        k = nterms;                % locations of model with all terms
    case 2
        modw = termsnotcontained ( fulltermlist, continuous, dohier );
        k = ( sum ( modw, 2 ) == nterms );
        TnotC = modw;
    case 3
        modw = true ( nterms );
        k = 1: nterms;
        TnotC = termsnotcontained ( fulltermlist, continuous, dohier );
end

% Gets a mdoel without the diagonal.
modwo = modw & ( eye ( nterms ) == 0 );


% STEP 7:  Fit each model, get its residual SS and d.f.
% for full model we already know the results
ssw ( k, : ) = repmat ( ssx, numel ( k ), 1 );
dfw ( k, : ) = repmat ( dfx, numel ( k ), 1 );

% Fits each model separately.
dfboth = cat ( 1, dfw, dfwo );

if doems
    Llist = cell ( nterms, ncomps );
end

% Considers interactions before their components for type-3 ss, so
% examines the terms in decreasing order of the number factors in the term.
if sstype == 3
    [ ~, sindices ] = sort ( sum ( fulltermlist, 2 ) );
    sindices = cat ( 1, sindices (:), sindices (:) );
else
    sindices = cat ( 2, 1: size ( termlist, 1 ), 1: size ( termlist, 1 ) )';
end

% Fits all required subsets of the full model.
fitsubmodels;

% STEP 7A:  Computes expected means squares, if requested.
if doems
    emsMat = makeemsmat ( dmat, cmat, nterms, Llist, termname );
end


% STEP 8:  Computes the sum of squares attributed to each term.
ssterm = max ( 0, ssw - sswo );
dfterm = dfw - dfwo;
ssterm ( dfterm == 0 ) = 0;


% STEP 9:  Computes the mean square for each term.
msterm = ssterm ./ max ( 1, dfterm );
mse    = sse .* ( dfe > 0 ) ./ max ( 1, dfe );

% STEP 9a:  Equates computed and expected mean squares, then solves for
%           variance component estimates
randomterm = find ( termlist * randomvar > 0 );
if isempty ( randomterm )
    msdenom = repmat ( mse, size ( msterm, 1 ), 1 );
    dfdenom = repmat ( dfe, size ( msterm, 1 ), 1 );
else
    [msdenom,dfdenom,varest,varci,txtdenom,txtems,denommat,rtnames] = ...
        getRandomInfo(msterm,dfterm,mse,dfe,emsMat,randomterm,...
        terminfo.tnames,alpha);
end


% STEP 10:  Computes the F statistic for each term and comparison.
fstat  = Inf ( size ( msterm ) );
t1     = ( msdenom > 0 );
t2     = ( msterm  > 0 );
fstat  ( t1 ) = msterm ( t1 ) ./ msdenom ( t1 );
fstat  ( ~( t1 | t2 ) ) = NaN;


% STEP 11:  Computes the p-value for each term and comparison.
pval   = NaN ( size ( fstat ) );
t      = ( dfdenom > 0 & dfterm > 0 );
pval   ( t ) = fpval ( fstat ( t ), dfterm ( t ), dfdenom ( t ) );
p      = pval;


% Reshapes the p-values to fit the original comparisons' shape.
pshape = mshape;
pshape (1) = size ( p, 1 );
p = reshape ( p, pshape );

% Marks the singular terms.
sing = bsxfun ( @lt, dfterm, terminfo.dfterm0 );
if ~isempty ( tnested )
    sing ( any ( tnested, 2 ), : ) = false;
end


% Creates the ANOVA table as a cell array, if requested.
if nargout >= 2 || isequal ( displayopt, 'on' )
    table = maketable;
end

% Returns additional information for multiple comparisons, if requested.
if nargout >= 3
    stats = makestats ( stats );  % nested function
end

% Returns the model terms, if requested.
if nargout > 3
    terms = termlist;
end



% ----- nested function
function fitsubmodels
   for j=length(sindices):-1:1
      % Find the next model index to fit
      k = sindices(j);
      
      % Look in unsorted arrays to see if we have already fit this model
      if j>nterms
         k0 = k+nterms;
      else
         k0 = k;
      end
      if dfboth(k0)~=-1
         continue
      end
      
      % Find the model with this index
      if (j > nterms)
         thismod = modwo(k, :);
      else
         thismod = modw(k, :);
      end
      
      % Get the design matrix for this model
      keepterms = find(thismod);
      clist = ismember(termname, [0 keepterms]);
      X = dmat2(:,clist);
      C = cmat(:,clist);

      % Fit this model
      [ ssx0, dfx0 ] = dofit ( X, data, C );
      
      % If this model is the "without" model for a term, compute L
      if doems && j>nterms
         if sstype==1
            prevterms = (1:nterms)<k;
         else
            prevterms = TnotC(k,:);
         end
         termscontained = ~prevterms;
         prevterms(k) = false;
         prevcols = ismember(termname, [0 find(prevterms)]);
         curcols = (termname == k);
         if sstype==3
            prevL = cat(1,Llist{termscontained});
         else
            prevL = [];
         end
         Llist{k} = getL(dmat, curcols, prevcols, prevL);
      end
      
      % Use these results for each term that requires them
      mod0 = repmat(thismod, nterms, 1);
      k = find(all(modw == mod0, 2));
      ssw(k,:) = repmat ( ssx0, numel ( k ), 1 );
      dfw(k,:) = repmat ( dfx0, numel ( k ), 1 );
      dfboth(k) = 0;
      k = find(all(modwo == mod0, 2));
      sswo(k,:) = repmat ( ssx0, numel ( k ), 1 );
      dfwo(k,:) = repmat ( dfx0, numel ( k ), 1 );
      dfboth(nterms+k) = 0;
   end
end

% ----- nested function
function stats = makestats(stats)
   stats.terms = termlist;
   stats.nlevels = varinfo.vlevels;
   stats.nlevels(continuous) = 1;
   stats.continuous = continuous;
   stats.vmeans = varinfo.vmeans;
   stats.termcols = [1; terminfo.termlength];
   stats.coeffnames = names;
   stats.vars = codes;
   stats.varnames = varinfo.varnames; % names of grouping variables
   stats.grpnames = varinfo.allnames; % names of levels of these variables
   if length(stats.resid)<length(nanrow)
      r(~nanrow) = stats.resid;
      r(nanrow) = NaN;
      stats.resid = r;
   end
   stats.vnested = vnested;

   % Extra info for models with random effects
   if doems
      stats.ems = emsMat;
      stats.denom = denommat;
      stats.dfdenom = dfdenom;
      stats.msdenom = msdenom;
      stats.varest = varest;
      stats.varci = varci;
      stats.txtdenom = txtdenom;
      stats.txtems = txtems;
      stats.rtnames = rtnames;
   else
      stats.ems = [];
      stats.denom = [];
      stats.dfdenom = [];
      stats.msdenom = [];
      stats.varest = [];
      stats.varci = [];
      stats.txtdenom = [];
      stats.txtems = [];
      stats.rtnames = [];
   end
end


% ----- nested function
function table = maketable
    
    % Initializes the table.
    table = zeros ( nterms + 2, 7, ncomps );
    
    % Writes the results for each term.
    table ( 2: nterms + 1, 2, : ) = ssterm;
    table ( 2: nterms + 1, 3, : ) = dfterm;
    table ( 2: nterms + 1, 4, : ) = sing;
    table ( 2: nterms + 1, 5, : ) = msterm;
    table ( 2: nterms + 1, 6, : ) = fstat;
    table ( 2: nterms + 1, 7, : ) = pval;
    
    % Writes the results for the error.
    table ( end, 2, : ) = sse;
    table ( end, 3, : ) = dfe;
    table ( end, 4, : ) = 0;
    table ( end, 5, : ) = mse;
    table ( end, 6, : ) = nan;
    table ( end, 7, : ) = nan;
    
    % Transforms the table to a cell.
    table = num2cell(table);
    table ( end, end - 1: end, : ) = { [] };
    
    % Writes the column headers.
    table ( 1, 1, : ) = cellstr ( getString ( message ( 'stats:anovan:ColHeadSource' ) ) );
    table ( 1, 2, : ) = cellstr ( getString ( message ( 'stats:anovan:ColHeadSumSq' ) ) );
    table ( 1, 3, : ) = cellstr ( getString ( message ( 'stats:anovan:ColHeadDF' ) ) );
    table ( 1, 4, : ) = cellstr ( getString ( message ( 'stats:anovan:ColHeadSingular' ) ) );
    table ( 1, 5, : ) = cellstr ( getString ( message ( 'stats:anovan:ColHeadMeanSq' ) ) );
    table ( 1, 6, : ) = cellstr ( getString ( message ( 'stats:anovan:ColHeadF' ) ) );
    table ( 1, 7, : ) = cellstr ( getString ( message ( 'stats:anovan:ColHeadProbGtF' ) ) );
    
    % Writes the row headers.
    for tindex = 1: nterms + 1
        table ( tindex + 1, 1, : ) = terminfo.tnames ( tindex );
    end
    
    % Add random effects information if required
    if ~isempty ( randomterm )
        
        % Adds eight new columns.
        ncols = size ( table, 2 );
        table ( :, end + 8, : ) = {[]};
        
        % Writes out the eight new headers.
        table ( 1, ncols + 1, : ) = cellstr ( getString ( message ( 'stats:anovan:RowHeadType' ) ) );
        table ( 1, ncols + 2, : ) = cellstr ( getString ( message ( 'stats:anovan:RowHeadExpectedMS' ) ) );
        table ( 1, ncols + 3, : ) = cellstr ( getString ( message ( 'stats:anovan:RowHeadMSDenom' ) ) );
        table ( 1, ncols + 4, : ) = cellstr ( getString ( message ( 'stats:anovan:RowHeadDfDenom' ) ) );
        table ( 1, ncols + 5, : ) = cellstr ( getString ( message ( 'stats:anovan:RowHeadDenomDefn' ) ) );
        table ( 1, ncols + 6, : ) = cellstr ( getString ( message ( 'stats:anovan:RowHeadVarEst' ) ) );
        table ( 1, ncols + 7, : ) = cellstr ( getString ( message ( 'stats:anovan:RowHeadVarLowerBnd' ) ) );
        table ( 1, ncols + 8, : ) = cellstr ( getString ( message ( 'stats:anovan:RowHeadVarUpperBnd' ) ) );
        
        % Marks the error as random.
        rterms = cat ( 1, randomterm, nterms + 1 );
        
        % Labels the random variables.
        table ( 2: nterms + 1, ncols + 1, : ) = cellstr ( getString ( message ( 'stats:anovan:RowHeadFixed' ) ) );
        table ( rterms + 1, ncols + 1, : )    = cellstr ( getString ( message ( 'stats:anovan:RowHeadRandom' ) ) );
        
        itxtems   = 1: size ( txtems,     1 );
        imsdenom  = 1: size ( msdenom,    1 );
        itxtdenom = 1: size ( txtdenom,   1 );
        
        table ( itxtems   + 1, ncols + 2, : ) = txtems;
        table ( imsdenom  + 1, ncols + 3, : ) = num2cell ( msdenom );
        table ( imsdenom  + 1, ncols + 4, : ) = num2cell ( dfdenom );
        table ( itxtdenom + 1, ncols + 5, : ) = txtdenom;
        table ( rterms    + 1, ncols + 6, : ) = num2cell ( varest );
        table ( rterms    + 1, ncols + 7, : ) = num2cell ( varci ( :, 1, : ) );
        table ( rterms    + 1, ncols + 8, : ) = num2cell ( varci ( :, 2, : ) );
    end
    
    % Adds a line with the total information.
    table ( end + 1, :, : ) = {[]};
    
    table ( end, 1, : ) = cellstr ( getString ( message ( 'stats:anovan:RowHeadTotal' ) ) );
    table ( end, 2, : ) = num2cell ( sst );
    table ( end, 3, : ) = num2cell ( ncases - 1 );
    table ( end, 4, : ) = num2cell (0);
    
    % Displays the results, if requested.
    if isequal ( displayopt, 'on' )
        
        % Doesn't show the results for more than one comparison.
        if ncomps > 1
            warning ( 'Too many comparisons. Disabling results'' table display.' )
            return
        end
        
        % Creates a reduced version of the table.
        dtable = table ( :, [ 1 2 3 5 6 7 ], : );
        
        % Marks the singular terms with a hash mark.
        [ sr, sc ] = find ( sing );
        dtable ( sr + 1, 1, sc ) = strcat ( { '# ' }, dtable ( sr + 1, 1, sc ) );
        
        % Sets the bottom-line message (type of sum-of-squares).
        switch sstype
            case 1
                cap = getString ( message ( 'stats:anovan:SequentialTypeISumsOfSquares' ) );
            case 2
                if dohier
                    cap = getString ( message ( 'stats:anovan:HierarchicalmodifiedTypeIISumsOfSquares' ) );
                else
                    cap = getString ( message ( 'stats:anovan:HierarchicalTypeIISumsOfSquares' ) );
                end
            otherwise
                cap = getString ( message ( 'stats:anovan:ConstrainedTypeIIISumsOfSquares' ) );
        end
        
        % Adds the legend to the bottom-line mssage.
        if any ( sing (:) )
            cap = [ cap '  ' getString(message('stats:anovan:TermsMarkedAreNotFullRank'))];
        end
        digits = [ -1 -1 -1 -1 2 4 ];
        
        % Sets the table title and header.
        title = getString ( message ( 'stats:anovan:NWayANOVA' ) );
        header = getString ( message ( 'stats:anovan:AnalysisOfVariance' ) );
        
        % Diaplays as many tables as comparisons.
        for cindex = 1: ncomps
            figh = statdisptable ( dtable ( :, :, cindex ), title, header, cap, digits, figure );
            set ( figh, 'HandleVisibility', 'callback' );
        end
    end
end % of nested function maketable

end % of main function


% --------------------------
function m = termsnotcontained(terms,continuous,dohier)
%TERMSNOTCONTAINED Creates a logical matrix indicating term containment
%   m(i,j)==1  iff  t(i) is not contained by t(j)

% For starters, omit continuous variables from this test
contnum = find(continuous);
t = terms;
t(:,contnum) = 0;

% The main test:  no overlap between terms included in t(i) and not in t(j)
m = (t*~t') > 0;

% Now consider continuous variables (X)
nterms = size(terms,1);
for j=1:length(contnum)
   v = terms(:,contnum(j));
   vv = repmat(v,1,nterms);
   if dohier
      % A cannot be contained in B if A has a higher power of X
      m = m | vv>vv';
   else
      % A cannot be contained in B if they don't match on X
      m = m | vv~=vv';
   end
end

% Set diagonals to 1 because we want proper containment
m(1:(nterms+1):end) = true;
end

% --------------------------
function v = makemodel(p,ng,vnested)
%MAKEMODEL Helper function to make model matrix
%    P = max model term size, NG = number of grouping variables
% or
%    P = vector of term codes
% VNESTED is a matrix specifying direct or indirect nesting.

% We want a matrix with one row per term.  Each row has a 1 for
% the variables participating in the term.  There will be rows
% summing up to 1, 2, ..., p.

if numel(p)==1
   % Create model matrix from a scalar max order value
   vgen = 1:ng;
   v = eye(ng);                      % linear terms
   for j=2:min(p,ng)
      c = nchoosek(vgen,j);          % generate column #'s with 1's
      nrows = size(c,1);             % generate row #'s
      r = repmat((1:nrows)',1,j);    %    and make it conform with c
      m = zeros(nrows,ng);           % create a matrix to get new rows
      m(r(:)+nrows*(c(:)-1)) = 1;    % fill in 1's
      v = cat ( 1, v, m );           % append rows
   end
else
   % Create model matrix from terms encoded as bit patterns
   nterms = length(p);
   v = zeros(nterms,ng);
   for j=1:nterms
      tm = p(j);
      while(tm)
         % Get last-numbered effect remaining
         lne = 1 + floor(log2(tm));
         tm = bitset(tm, lne, 0);
         v(j,lne) = 1;
      end
   end
end

% Remove terms forbidden by nesting
[nestee,nester] = find(vnested);
bad = any(v(:,nestee) & v(:,nester), 2);
v(bad,:) = [];

end

% --------------------------
function [ssx,dfx,dmat2,y2,stats] = dofit(dmat,y,cmat,sst,mu)
%DOFIT Do constrained least squares fit and reduce data for subsequent fits

% Gets the number of comparisons.
ncomps = size ( y, 2 );

% Find the null space of the constraints matrix
[Qc,Rc,~] = qr(cmat');
pc = Rrank(Rc);
Qc0 = Qc(:,pc+1:end);

% Do qr decomposition on design matrix projected to null space
Dproj = dmat*Qc0;
[Qd,Rd,Ed] = qr(Dproj,0);
clear Dproj

dfx = Rrank(Rd);
Qd = Qd(:,1:dfx);
Rd = Rd(1:dfx,1:dfx);

% Fit y to design matrix in null space
y2 = Qd' * y;            % rotate y into that space
zq = Rd \ y2;            % fit rotated y to projected design matrix

z = zeros(length(Ed),ncomps); % coefficient vector extend to full size ...
z(Ed(1:dfx),:) = zq;       % ... and in correct order for null space

b = Qc0 * z;             % coefficients back in full space
% ssx = norm(y2)^2;        % sum of squares explained by fit
ssx = sqrt ( sum ( y2 .^ 2, 1 ) ) .^ 2;

% Return reduced design matrix if requested
if nargout > 2
   dmat2 = Qd' * dmat;
end

% Prepare for multiple comparisons if requested
if nargout >= 3
    stats.source = 'anovan';
    
    % Get residuals
    % yhat = Dproj * z;
    yhat = Qd * Rd * zq;
    stats.resid = y - yhat;
    
    % Calculate coefficients, then adjust for previously removed mean
    t = b;
    if ~isempty ( t )
        
        % Undoes the mean adjustment.
        t ( 1, : ) = t ( 1, : ) + mu;
    end
    
    % Stores the data.
    stats.coeffs = t;
    
    RR = zeros(dfx,length(Ed));
    RR(:,Ed(1:dfx)) = Rd;
    stats.Rtr = RR';
    [~,rr,ee] = qr(dmat,0);
    pp = Rrank (rr );
    rowbasis = zeros(pp,size(rr,2));
    rowbasis(:,ee) = rr(1:pp,:);
    stats.rowbasis = rowbasis;
    stats.dfe = size ( y, 1 ) - dfx;
    stats.mse = ( sst - ssx ) / max ( 1, stats.dfe );
    stats.nullproject = Qc0;
end

% Sets the degrees of freedom in matrix form.
dfx = repmat ( dfx, 1, ncomps );
end

% ---------------------
function p = Rrank(R,tol)
%RRANK Compute rank of R factor after qr decomposition
if (min(size(R))==1)
   d = abs(R(1,1));
else
   d = abs(diag(R));
end
if nargin<2
   tol = 100 * eps(class(d)) * max(size(R));
end
p = sum(d > tol*max(d));
end

% ----------------------
function ems = getems(L,R,termnums)
%GETEMS Final step in expected mean square computation
%   The L array is an array of constraints on the coefficients,
%   R is the R factor in a QR decomposition X, and termnums is a
%   vector mapping R column numbers to term numbers.  The output is
%   a vector of expected mean square coefficients with one element
%   per term.  The first element is for the constant term, and
%   the last is for the error term.  See "Computing Expected Mean
%   Squares" by Goodnight and Speed, 1978.

% In G&S notation, form the matrix L*pinv(X'*X)*L'.
dfL = size(L,1);
temp = L/R;
LML = temp*temp';

% Create the Cholesky decomposition U and then compute U\L.
U = cholcov(LML);
C = U'\L;

% Sum the squared values of the columns other than the column for
% the constant term, and assign the sums to terms.
sumsq = sum(C(:,2:end).^2,1);
cols = termnums(2:end);
rows = ones(size(cols));
ems = accumarray([rows(:),cols(:)], sumsq) / dfL;
ems(:,end+1) = 1;  % add entry for error term
end

% -----------------------------
function L = getL(X, curcols, incols, L2)
%GETL Get L hypothesis matrix for the current term after adjusting for
%     terms that remain "in" for this test, making the whole thing
%     orthogonal to the matrix L2 containing hypotheses for terms that
%     contain the current term.
%
%        Type   incols                  L2
%        ----   --------------------    ----------------
%        1      lower-numbered terms    empty
%        2      non-containing terms    empty
%        3      non-containing terms    containing terms

x1 = X(:,curcols);     % the current term
x0 = X(:,incols);      % terms that remain "in" when testing x1

% Find hypothesis for the current term.  First remove effects of terms it does
% not contain (these terms remain in the model when testing the current term).
[Q0,R,~] = qr(x0,0);
tol = 100 * max(size(x0)) * eps(R(1,1));
if min(size(R))>1
   p = sum(abs(diag(R)) > tol);
else
   p = 1;
end

Q0 = Q0(:,1:p);
adjx1 = x1 - Q0 * (Q0' * x1);

% Now find a full-rank basis for the adjusted effect.
[Q2,R,~] = qr(adjx1,0);
tol = 100 * max(size(R)) * eps(R(1,1));
if min(size(R))>1
   p = sum(abs(diag(R)) > tol);
else
   p = 1;
end
Q2 = Q2(:,1:p);
L = Q2' * X;

% Make L rows orthogonal to the L rows for containing terms.
if ~isempty(L2)
   L2tr = L2';
   [Q,~,~] = qr(L2tr,0);
   L = (L' - Q * (Q' * L'))';

   % Make L full rank
   [~,R,~] = qr(L,0);
   tol = 100 * max(size(L)) * eps(R(1,1));
   if min(size(R))>1
      p = sum(abs(diag(R)) > tol);
   else
      p = 1;
   end
   if p<size(L,1)
      L = L(1:p,:);
   end
end
end

% -----------------------------
function [msDenom,dfDenom,varest,varci,txtdenom,txtems,denommat,rtnames] = ...
                      getRandomInfo(msTerm,dfTerm,msErr,dfErr,ems,...
                                    randomterms,tnames,alpha)
%GETRANDOMINFO Get info for random effects, such as denominator for F tests

nterms   = size ( msTerm, 1 );
ncomps   = size ( msTerm, 2 );
msDenom  = nan  ( nterms, ncomps );
dfDenom  = nan  ( nterms, ncomps );
txtdenom = cell ( nterms, ncomps );
txtems   = cell ( nterms + 1, ncomps );
txtems ( end, : ) = cellstr ( getString ( message ( 'stats:anovan:VError' ) ) );

% For convenience, combine error info with term info
msTerm = cat ( 1, msTerm, msErr );
dfTerm = cat ( 1, dfTerm, dfErr );
randomterms = cat ( 1, randomterms (:), size ( msTerm, 1 ) );

% Find the matrix giving expected mean squares for random rows, and
% compute the coefficients of these expected values.
randrows = ems ( randomterms, : );
emsmat = ems ( 1: end - 1, : );
emsmat ( 1: size ( emsmat, 1 ) + 1: end ) = 0;
[ q, r, e ] = qr ( randrows', 0 );
p = Rrank ( r );
B = r ( 1: p, : ) \ ( q ( :, 1: p )' * emsmat' );

% Loop over all terms
for j = 1: nterms
   % First compute F statistic denominator information
   % See which terms are involved in this combination
   b = B ( :, j );
   tol = sqrt ( eps ( class ( b ) ) ) * max ( abs ( b ) );
   nonzeroterms = find ( abs ( b ) >= tol );
   indices = e ( nonzeroterms );
   tnums = randomterms ( indices );
   
   % If it's a single term, uses it exactly
   if length ( nonzeroterms ) == 1
      msDenom  ( j, : ) = msTerm ( tnums, : );
      dfDenom  ( j, : ) = dfTerm ( tnums, : );
      txtdenom ( j, : ) = cellstr ( getString ( message ( 'stats:anovan:MSParen', tnames { tnums } ) ) );
      
   % Otherwise uses linear combination and Satterthwaite's approximation.
   else
      linprod = bsxfun ( @times, b ( nonzeroterms ), msTerm ( tnums, : ) );
      msDenom ( j, : ) = sum ( linprod, 1 );
      dfDenom ( j, : ) = msDenom ( j, : ) .^ 2 ./ sum ( ( linprod .^ 2 ) ./ dfTerm ( tnums, : ) );
      dfDenom ( j, any ( dfTerm ( tnums, : ) == 0 ) ) = nan;
      txtdenom ( j, : ) = cellstr ( linprodtxt ( b ( nonzeroterms ), tnames ( tnums ), getString ( message ( 'stats:anovan:MS' ) ) ) );
   end
   
   % Next get an expression for the expected mean square
   prefix = repmat ( getString ( message ( 'stats:anovan:Q' ) ), nterms + 1, 1 ); % quadratic form for fixed effects
   prefix ( randomterms, 1 ) = getString ( message ( 'stats:anovan:V' ) );     % variance for random effects
   txtems { j } = linprodtxt ( ems ( j, : ), tnames, prefix );
end

% Get the variance component estimates
[ varest, varci ] = varcompest ( ems, randomterms, msTerm, dfTerm, alpha );

% Package up the results
denommat = B';
rtnames = tnames ( randomterms );
end

% -----------------------------------------
function [varest,varci] = varcompest(ems,Rterms,msTerm,dfTerm,alpha)
%VARCOMPEST Estimate variance components from obs. & exp. mean squares

% Gets the number of therms excluding the error.
nterms = size ( msTerm, 1 ) - 1;
ncomps = size ( msTerm, 2 );

Fterms = find ( ~ismember ( 1: nterms, Rterms ) );
Fterms = Fterms (:);

% Get A and B so that expected means square is Av+Bq
% where v is the variances and q is quadratic forms in the fixed effects
B = ems ( Rterms, Fterms );
A = ems ( Rterms, Rterms );
const = msTerm ( Rterms, : ) - B * msTerm ( Fterms, : );
Ainv = pinv ( A );          % should be well-conditioned and easily invertible
varest = Ainv * const;

% Get confidence bounds for the variance components
varci = nan ( size ( varest, 1 ), 2, ncomps );
if all ( dfTerm (:) > 0 )
   L = msTerm .* dfTerm ./ chi2inv ( 1 - alpha / 2, dfTerm );
   U = msTerm .* dfTerm ./ chi2inv ( alpha / 2, dfTerm );
   AinvB = Ainv * B;
   
   for j = find ( any ( varest > 0, 2 )' )
      Rcoeffs = Ainv  ( j, : )';
      Fcoeffs = AinvB ( j, : )';
      t = ( Rcoeffs > 0 );
      s = ( Fcoeffs > 0 );
      cL = zeros ( 1, ncomps );
      cU = zeros ( 1, ncomps );
      if any ( t )
         cL = cL + sum ( bsxfun ( @times, Rcoeffs ( t ),  L ( Rterms ( t ), : ) ), 1 );
         cU = cU + sum ( bsxfun ( @times, Rcoeffs ( t ),  U ( Rterms ( t ), : ) ), 1 );
      end
      if any ( ~t )
         cL = cL + sum ( bsxfun ( @times, Rcoeffs ( ~t ), U ( Rterms ( ~t ), : ) ), 1 );
         cU = cU + sum ( bsxfun ( @times, Rcoeffs ( ~t ), L ( Rterms ( ~t ), : ) ), 1 );
      end
      if any ( s )
         cL = cL + sum ( bsxfun ( @times, Fcoeffs ( s ),  L ( Fterms ( s ), : ) ), 1 );
         cU = cU + sum ( bsxfun ( @times, Fcoeffs ( s ),  U ( Fterms ( s ), : ) ), 1 );
      end
      if any ( ~s )
         cL = cL + sum ( bsxfun ( @times, Fcoeffs ( ~s ), U ( Fterms ( ~s ), : ) ), 1 );
         cU = cU + sum ( bsxfun ( @times, Fcoeffs ( ~s ), L ( Fterms ( ~s ), : ) ), 1 );
      end
      varci ( j, 1, : ) = max ( 0, cL );
      varci ( j, 2, : ) = max ( 0, cU );
   end   
end
end

% ------------------------------------------------
function txt=linprodtxt(coeffs,names,prefix)
%LINPRODTXT Create a text representing a linear combination of some things

txt = '';
plustxt = '';

% If things will display as 0 or 1, make them exact
tol = sqrt(eps(class(coeffs))) * max(abs(coeffs));
coeffs(abs(coeffs)<tol) = 0;
coeffs(abs(coeffs-1)<tol) = 1;
coeffs(abs(coeffs+1)<tol) = -1;

% Add each component of the linear combination
preftxt = prefix;
for i=1:length(coeffs)
   bi = coeffs(i);
   if bi~=0
      % Get the sign
      if bi<0
         signtxt = '-';
      else
         signtxt = plustxt;
      end
      
      % Get the coefficient and multiplication symbol, unless it's one
      if abs(bi)==1
         coefftxt = '';
      else
         coefftxt = sprintf('%g*',abs(bi));
      end
      
      % If each row has its own prefix, get this one
      if size(prefix,1)>1
         preftxt = deblank(prefix(i,:));
      end
      
      txt = sprintf('%s%s%s%s(%s)',txt,signtxt,coefftxt,preftxt,names{i});
      plustxt = '+';
   end
end
end

% -------------------
function [varnames,termlist,randomvar,sstype,continuous,vnested] = ...
        argcheck(ng,varnames,model,randomvar,sstype,alpha,continuous,vnested)
% variable names
if isempty(varnames)
   varnames = cellstr([repmat('X',ng,1) strjust(num2str((1:ng)'),'left')]);
else
   if (iscell(varnames)), varnames = varnames(:); end
   if (size(varnames,1) ~= ng)
      error(message('stats:anovan:VarNamesSizeMismatch', ng));
   elseif (~(ischar(varnames) || iscellstr(varnames)))
      error(message('stats:anovan:BadVarNames'));
   elseif (ischar(varnames))
      varnames = cellstr(varnames);
   end
end

% sum-of-squares type
if isempty(sstype)
   sstype = 3;
elseif ~isscalar(sstype) || (sstype ~= 1 && sstype ~= 2 && sstype ~=  3 && lower(sstype) ~= 'h')
   error(message('stats:anovan:BadSumSquares'));
end

% continuous must be a logical vector indicating continuous variables
if islogical(continuous) && length(continuous)==ng && numel(continuous)==ng
   continuous = continuous(:);
elseif isnumeric(continuous) && all(ismember(continuous(:),1:ng))
   continuous = ismember(1:ng,continuous);
else
   error(message('stats:anovan:BadContinuous', ng));
end

% nested must be a binary matrix
if ~isempty(vnested)
    if ~isequal(size(vnested),[ng ng]) || ~all(vnested(:)==0 | vnested(:)==1)
        error(message('stats:anovan:BadNested', ng, ng));
    end
    
    % find all indirect nesting
    vnested = nest_d2i(vnested);
end

% alpha (confidence level is 100*(1-alpha)%
if ~isnumeric(alpha) || numel(alpha)~=1 || alpha<=0 || alpha>=1
   error(message('stats:anovan:BadAlpha'));
end

% Randomvar may be a logical vector indicating the random variables
if islogical(randomvar) && length(randomvar)==ng && numel(randomvar)==ng
   randomvar = randomvar(:);
elseif isnumeric(randomvar) && all(ismember(randomvar(:),1:ng))
   % or it may be an index vector -- convert to logical
   randomvar = ismember((1:ng)', randomvar(:));
else
   error(message('stats:anovan:BadRandomVar', ng));
end
if any(randomvar(continuous))
   error(message('stats:anovan:ContinuousRandomVar'));
end

% model
if isempty(model)
   model = 1;
elseif ~isnumeric(model)
   if strcmp(model, 'linear') || strcmp(model, 'l')
      model = 1;
   elseif strcmp(model, 'interaction') || strcmp(model, 'i')
      model = 2;
   elseif strcmp(model, 'full') || strcmp(model, 'f')
      model = ng;
   else
      error(message('stats:anovan:UnrecognizedModel'));
   end
elseif min(size(model))>1 || ...
       (isequal(size(model),[1 ng]) && all(ismember(model,0:1)))
   % Matrix form of input
   if (any(model~=floor(model)|model<0))
      error(message('stats:anovan:BadModelMatrix'));
   elseif size(model,2)~=ng
      error(message('stats:anovan:ModelMatrixColumnSizeMismatch', ng));
   end
else
   if (any(model~=floor(model)|model<=0))
      % Original vector form, no longer advertised
      error(message('stats:anovan:BadModelVector'));
   end
   model = model(:);
end

% convert model to numeric matrix list of terms
if (length(model)==1) || (any(size(model)==1) && size(model,2)~=ng)
   % Convert from scalar or coded integer to matrix form
   termlist = makemodel(model, ng, vnested);
else
   termlist = model;
end

% look for invalid model terms
t = all(termlist==0,2);
termlist(t,:) = [];      % remove constant term -- it's implicit
if any(any(termlist(:,~continuous)>1))
    error(message('stats:anovan:ModelNotLinear'));
end
if ~isempty(vnested)
   % Must not have A(B) and B(A)
   if any(any(vnested & vnested'))
       error(message('stats:anovan:BadNesting'));
   end
   
   % Must not have an interaction between B and A(B)
   [nestee,nester] = find(vnested);
   bad = any(termlist(:,nestee) & termlist(:,nester), 2);
   if any(bad)
      error(message('stats:anovan:BadModelNestedInteraction'));
   end
end
end

% ------------------
function  emsMat = makeemsmat(dmat, cmat, nterms, Llist, termname)

% Factor the design matrix X appended with the constraints matrix C,
% getting a full-rank R factor.
xc = [dmat; cmat];
[~,Rxc,Exc] = qr(xc,0);
tol = sqrt(eps(class(Rxc)));
rankRxc = Rrank(Rxc,tol);
Rxc0 = zeros(rankRxc,size(Rxc,2));
Rxc0(:,Exc) = Rxc(1:rankRxc,:);

% Compute expected mean squares for each term.
emsMat = zeros(nterms,0);
for k=1:nterms
   termems = getems(Llist{k},Rxc0,termname);
   emsMat(k,1:length(termems)) = termems;
end
emsMat(nterms+1,end) = 1;       % for error term
end

% ----------------------
function [codes,names] = convertstats(stats,terminfo,varinfo,continuous)

allnames = varinfo.allnames;
varnames = varinfo.varnames;
ng = length ( varnames );
fullcoef = stats.coeffs;
codes = zeros ( size ( fullcoef, 1 ), ng );
names = cell ( size ( fullcoef, 1 ), 1 );
base = 1;
names {1} = getString ( message ( 'stats:anovan:Constant' ) );

% compute names
for j=1:length(terminfo.levelcodes)
   M = terminfo.levelcodes{j};
   varlist = terminfo.termvars{j};
   v1 = varlist(1);
   vals1 = allnames{v1};
   nrows = size(M,1);
   codes(base+(1:nrows),varlist) = M;
   for k=1:nrows
      if continuous(v1)
          nm = varnames{v1};
      else
          nm = sprintf('%s=%s',varnames{v1},vals1{M(k,1)});
      end
      for i=2:length(varlist)
         v2 = varlist(i);
         vals2 = allnames{v2};
         if continuous(v2)
             nm2 = varnames{v2};
         else
             nm2 = sprintf('%s=%s',varnames{v2},vals2{M(k,i)});
         end
         nm = sprintf('%s * %s',nm,nm2);
      end
      names{base+k} = nm;
   end
   base = base+nrows;
end
end

% --------------------------------
function terminfo = makedummyterms(sindex,termlist,varinfo,continuous,tnested)

ncols = 1;
[nterms,nfactors] = size(termlist);

termdum = cell(nterms, 1);      % dummy variable = design matrix cols
termconstr = cell(nterms,1);    % constraints to make each term well defined
levelcodes = cell(nterms, 1);   % codes for levels of each M row
tnames = cell(nterms, 1);       % name for entire term, e.g. A*B
dfterm0 = zeros(nterms, 1);     % nominal d.f. for each term
termvars = cell(nterms, 1);     % list of vars in each term
termlength = zeros(size(sindex));% length of each term (number of columns)

auxtermlist = zeros(0,nfactors);  % may need to create temp terms for nesting
auxtermdum = cell(0,1);           % may need their dummy variables

isnested = false;

% For each term,
for j=1:nterms
   % Get dummy columns, var list, term name, etc
   sj = sindex(j);
   tm = termlist(sj,:);
   if ~isempty(tnested)
      isnested = any(tnested(sj,:));
   end
   [tdum,vars,tn,df0,tconstr] = maketerm(tm,isnested,varinfo,j,sindex,...
                                         termlist,termdum,termconstr,tnames,dfterm0);

   % Store this term's information
   k = size(tdum, 2);
   termlength(sindex(j),1) = k;
   ncols = ncols + k;
   termdum{sj} = tdum;
   termvars{sj} = vars;
   levelcodes{sj} = fliplr(fullfact(varinfo.vlevels(vars(end:-1:1))));
   tnames{sj,1} = tn;
   dfterm0(sj) = df0;

   % For a nested term, figure out the constraint now
   if isnested
      % First get dummy variables for this term and its nesters
      if any(tm(continuous))
         % Ignore any continuous variable contributions to this term
         tmcat = tm;                % this term, cat predictors only
         tmcat(continuous) = 0;
         ntmcat = termlist(tnested(sj,:),:);
         ntmcat(:,continuous) = 0;  % nesting terms, cat predictors only
         [Ydum,auxtermlist,auxtermdum] =findtermdum(tmcat,termlist,termdum,...
                              auxtermlist,auxtermdum,varinfo);
         Ydum = Ydum{1};
         [Xdum,auxtermlist,auxtermdum] =findtermdum(ntmcat,termlist,termdum,...
                              auxtermlist,auxtermdum,varinfo);
      else
         % No continuous variables, so work with the whole term
         Ydum = tdum;
         Xdum = termdum(tnested(sj,:)>0);
      end

      % Constrain the part of this term in its nesters to be 0
      tconstr = gettermconstr(Ydum, Xdum);
   end
   termconstr{sj} = tconstr;
end
tnames{length(tnames)+1,1} = getString(message('stats:anovan:Error'));

% Package up term info into a structure
terminfo.termdum = termdum;
terminfo.termconstr = termconstr;
terminfo.levelcodes = levelcodes;
terminfo.tnames = tnames;
terminfo.dfterm0 = dfterm0;
terminfo.termvars = termvars;
terminfo.termlength = termlength;
end

% -----------------------------------
function [y,allgrps,nanrow,varinfo] = removenans(y,group,continuous)

% Find NaNs among response and group arrays
ncases = size ( y, 1 );
nanrow = any ( isnan ( y (:, : ) ), 2 );
ng = length(group);
for j=1:ng
   gj = group{j};
   if (size(gj,1) == 1), gj = gj(:); end
   if (size(gj,1) ~= ncases)
      error(message('stats:anovan:GroupVarSizeMismatch', j, ncases));
   end
   if (ischar(gj)), gj = cellstr(gj); end
   if ~isvector(gj)
       error(message('stats:anovan:GroupNotVectorOrCharArray'));
   end

   group{j} = gj;
   if (isnumeric(gj))
      nanrow = (nanrow | isnan(gj));
   elseif isa(gj,'categorical')
      nanrow = (nanrow | isundefined(gj));
   else
      nanrow = (nanrow | strcmp(gj,''));
   end
end

% Remove rows with NaN anywhere
y ( nanrow, : ) = [];
ncases = size ( y, 1 );
ng = length(group);
dfvar = zeros(ng,1);
allgrps = zeros(ncases, ng);
allnames = cell(ng,1);

% Get arrays describing the groups
for j=1:ng
   gj = group{j};
   gj(nanrow,:) = [];
   group{j} = gj;
   if continuous(j)
      dfvar(j) = 1;
      allgrps(:,j) = gj;
      allnames{j} = {''};
   else
      [gij,gnj] = grp2idx(gj);
      nlevels = size(gnj,1);
      dfvar(j) = nlevels - 1;
      allgrps(:,j) = gij;
      allnames{j} = gnj;
   end
end

% The df and allnames information does not yet reflect nesting.  We will
% fix up allnames for nested variables later, and not use df for them.

varinfo.df = dfvar;
varinfo.allnames = allnames;
end

% -----------------------------
function varinfo = makedummyvars(continuous,allgrps,vnested,varinfo)
% Create main effect term values for each variable

ng = size(allgrps,2);
vdum = cell(ng,1);      % dummy terms for variables
vconstr = cell(ng,1);   % constraints
vlevels = ones(ng,1);   % levels corresponding to each dummy term
vmeans = zeros(ng,1);   % variable means, used for continuous factors
allnames = varinfo.allnames;

for j=1:ng
   if continuous(j)
      % For continuous variables, use the variable itself w/o constraints
      vdum{j} = allgrps(:,j);
      vconstr{j} = ones(0,1);
      vmeans(j) = mean(allgrps(:,j));
   elseif isempty(vnested) || ~any(vnested(j,:))
      % Create dummy variable arrays for each grouping variable
      % using a sum-to-zero constraint
      vdum{j} = idummy(allgrps(:,j), 3);
      nlevgrp  = size(vdum{j},2);
      vconstr{j} = ones(1,nlevgrp);  % sum-to-zero constraint for categorical
      vlevels(j) = nlevgrp;
   else
      % Create dummy variable arrays for a nested variable
      nesternums = find(vnested(j,:));
      allvars = allgrps(:,[nesternums, j]);
      [ugrps,~,grpnum] = unique(allvars,'rows');
      vdum{j} = idummy(grpnum,3);
      vconstr{j} = []; % to be computed later
      vlevels(j) = size(vdum{j},2);
      
      % Fix up names for nested variables
      levelnames = varinfo.allnames{j};
      namesj = cell(vlevels(j),1);
      nester1names = allnames{nesternums(1)};
      for rnum = 1:vlevels(j)
          nesterlist = nester1names{ugrps(rnum,1)};
          for k = 2:length(nesternums)
              nesterknames = allnames{nesternums(k)};
              nesterlist = sprintf('%s,%s',nesterlist,...
                                   nesterknames{ugrps(rnum,k)});
          end
          namesj{rnum} = sprintf('%s(%s)',levelnames{ugrps(rnum,end)},...
                                          nesterlist);
      end
      varinfo.allnames{j} = namesj;
   end
end

varinfo.vdum = vdum;
varinfo.vconstr = vconstr;
varinfo.vlevels = vlevels;
varinfo.vmeans = vmeans;
end

% -----------------------------
function [dmat,cmat,termname] = makedesignmat(terminfo,n)

termlength = terminfo.termlength;
ncols = sum(termlength);

nconstr = sum(cellfun('size',terminfo.termconstr,1));
dmat = ones(n, ncols+1);      % to hold design matrix
cmat = zeros(nconstr,ncols);  % to hold constraints matrix
cbase = 0;                    % base from which to fill in cmat
termname = zeros(ncols,1);
termstart = cumsum([2; termlength(1:end-1)]);
termend = termstart + termlength - 1;
for j=1:length(termlength)
   clist = termstart(j):termend(j);
   dmat(:, clist) = terminfo.termdum{j};
   C = terminfo.termconstr{j};
   nC = size(C,1);
   cmat(cbase+1:cbase+nC,clist) = C;
   termname(clist) = j;
   cbase = cbase + nC;
end
end


% -----------------------------
function tnested = gettermnesting(fulltermlist,vnested,continuous)
% Create matrix with (i,j) indicating if term i is nested in term j

% No nested variables implies no nested terms
if isempty(vnested)
    tnested = [];
    return
end

% Work with categorical and continuous terms separately
ctermlist = fulltermlist(:,continuous);
fulltermlist(:,continuous) = 0;

nterms = size(fulltermlist,1);
tnested = zeros(nterms);

% If A(B), then any term containing A may be nested in one containing B
[nestee,nester] = find(vnested);
for j=1:length(nester)
    hasnestee = (fulltermlist(:,nestee(j))>0);
    hasnester = (fulltermlist(:,nester(j))>0);
    tnested(hasnestee,hasnester) = 1;
end

% But terms are not nested within themselves
tnested(1:nterms+1:end) = 0;

% And there is no nesting relationship if the continuous vars don't match
for j=1:size(ctermlist,2)
    c = repmat(ctermlist(:,j),1,nterms);
    tnested = tnested & (c == c');
end

% And the full representation for a nesting term must be a subset of the
% full representation for a nested term
[nestee,nester] = find(tnested);
for j=1:length(nester)
    if any(fulltermlist(nestee(j),:) < fulltermlist(nester(j),:))
        tnested(nestee(j),nester(j)) = 0;
    end
end

end

% ----------------------------
function constr = gettermconstr(tdum,nesterdum)

% Create X matrix containing design matrix columns of nester terms
if iscell(nesterdum)
    Xnester = cat(2,nesterdum{:});
else
    Xnester = nesterdum;
end
nnested = size(Xnester,2);

% Find unique combinations of nested and nesting terms
Xboth = [Xnester tdum];
Uboth = unique(Xboth,'rows');

% Find constraints to force the part of the nested terms in the column
% space of the nesting terms to be zero
[Q,R,~] = qr(Uboth(:,1:nnested));
rnk = Rrank(R);
Q = Q(:,1:rnk);
constr = Q*(Q'*Uboth(:,nnested+1:end));

% Remove redundant constraints
[~,R,E] = qr(constr',0);
rnk = Rrank(R);
constr = constr(E(1:rnk),:);
end


% ----------------------------
function nestednames = makenestednames(varnames,vnested)

nestednames = varnames;

nestedvars = find(any(vnested,2));
for j=1:length(nestedvars)
    k = nestedvars(j);
    thelist = sprintf('%s,',varnames{vnested(k,:)>0});
    nestednames{k} = sprintf('%s(%s)',varnames{k},thelist(1:end-1));
end
end

% ----------------------------
function fulltermlist = showtermnesting(termlist,vnested)

fulltermlist = termlist;

[nestee,nester] = find(vnested);
for j=1:length(nester)
    t = fulltermlist(:,nestee(j))>0;
    fulltermlist(t,nester(j)) = 1;
end
end


% ----------------------------
function [tdum,vars,tn,df0,tconstr] = maketerm(tm,isnested,varinfo,j,sindex,termlist,...
                                             termdum,termconstr,tnames,dfterm0)
% Make term info such as dummy vars, name, constraints

% Loop over elements of the term
df0 = 1;
tdum = [];         % empty term so far
tconstr = 1;       % empty constraints so far
tn = '';           % blank name so far
vars = find(tm);   % list of variables making up terms
pwrs = tm(vars);   % and their exponents
tm = (tm>0);       % and a boolean mask for them
for varidx = 1:length(vars)
   % Process each variable participating in this term
   varnum = vars(varidx);          % variable name
   thispwr = pwrs(1);              % power of this variable
   tm(varnum) = 0;                 % term without this variable
   pwrs(1) = [];                   % powers without this variable
   df0 = df0 * varinfo.df(varnum); % d.f. so far

   % Combine its dummy variable with the part computed so far
   G = varinfo.vdum{varnum};       % dummy vars for this grouping var
   thisname = varinfo.varnames{varnum};
   if thispwr>1
       G = G .^ thispwr;
       thisname = sprintf('%s^%d',thisname,thispwr);
   end
   tdum = termcross(G,tdum);   % combine G into term dummy vars

   % Construct the term name and constraints matrix
   if nargout>1
      if (isempty(tn))
         tn = thisname;
         if ~isnested
            tconstr = varinfo.vconstr{varnum};
         end
      else
         tn = strcat ( tn, '*', thisname );
         if ~isnested
            tconstr = [kron(varinfo.vconstr{varnum},eye(size(tconstr,2)));
                       kron(eye(length(varinfo.vconstr{varnum})),tconstr)];
         end
      end
   end

   if varidx<length(vars) && j>1
      % If the rest of this term is computed, take advantage of that
      prevterms = termlist(sindex(1:j-1),:);
      oldtm = find(all(prevterms(:,tm) == repmat(pwrs,size(prevterms,1),1),2));
      oldtm = oldtm(~any(prevterms(oldtm,~tm),2));
      if ~isempty(oldtm)
         k = sindex(oldtm(1));
         tdum = termcross(termdum{k}, tdum);
         if nargout>1
            if ~isnested
               oconstr = termconstr{k};
               tconstr = [kron(tconstr,              eye(size(oconstr,2)));
                          kron(eye(size(tconstr,2)), oconstr)];
            end
            tn = strcat ( tn, '*', tnames {k} );
            df0 = df0 * dfterm0(k);
         end
         break;
      end
   end
end
if (isempty(tn)), tn = getString(message('stats:anovan:Constant')); end
end

% ----------------------------
function [dum,auxtermlist,auxtermdum] =findtermdum(tms,termlist,termdum,...
                   auxtermlist,auxtermdum,varinfo)

 
nterms = size(tms,1);
dum = cell(nterms,1);
for j=1:nterms
    tm = tms(j,:);

    % Try to find this term among the terms in the model
    oldtm = find(all(termlist == repmat(tm,size(termlist,1),1),2));
    if ~isempty(oldtm)
        dum{j} = termdum{oldtm(1)};
        continue
    end
        
    % Try to find this term among the auxiliary terms
    oldtm = find(all(auxtermlist == repmat(tm,size(auxtermlist,1),1),2));
    if ~isempty(oldtm)
        dum{j} = auxtermdum{oldtm(1)};
        continue
    end
    
    % Create an auxiliary term
    dum{j} = maketerm(tm,[],varinfo,0);
end

% Apends the new data.
auxtermdum  = cat ( 1, auxtermdum,  dum );
auxtermlist = cat ( 1, auxtermlist, tms );
end

% ---------------------------
function m = nest_d2i(nested)
% Take matrix specifying at least direct nesting relationships, and
% return matrix that includes all indirect nesting relationships
m = nested;
n = size(m,1);
m0 = zeros(n);
while(~isequal(m,m0))
   m0 = m;
   m = double(m | ((m*m) > 0));
end
end

function p = fpval(x,df1,df2)
%FPVAL F distribution p-value function.
%   P = FPVAL(X,V1,V2) returns the upper tail of the F cumulative distribution
%   function with V1 and V2 degrees of freedom at the values in X.  If X is
%   the observed value of an F test statistic, then P is its p-value.
%
%   The size of P is the common size of the input arguments.  A scalar input  
%   functions as a constant matrix of the same size as the other inputs.    
%
%   See also FCDF, FINV.

%   References:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.6.

%   Copyright 2010 The MathWorks, Inc. 
%   $Revision: 1.1.8.3 $  $Date: 2010/11/08 02:37:46 $

if nargin < 3, 
    error(message('stats:fpval:TooFewInputs')); 
end

xunder = 1./max(0,x);
xunder(isnan(x)) = NaN;
p = fcdf(xunder,df2,df1);

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
end

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
d = zeros ( n, ncols );

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
end