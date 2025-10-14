function pval = my_stats ( data, group, covariates, anova, parametric, nperms )


% Initializes the input variables.
if nargin < 6, nperms     = 0;  end
if nargin < 5, parametric = 1;  end
if nargin < 4, anova      = 0;  end
if nargin < 3, covariates = []; end

% Gets the number of grouping variables and of covariates, if any.
ngroup = size ( group, 2 );
ncov   = size ( covariates, 2 );

% If there are covariates the test must be ANOVA.
if ncov && ~parametric
    error ( 'Only parametric ANOVA can include covariates.' );
end

if ncov && ~anova
    warning ( 'Forcing the test to ANOVA to include covariates.' );
    anova = 1;
end

% If more than one grouping variable allows only ANOVA.
if ngroup > 1 && ~parametric
    error ( 'Only parametric ANOVA can use more than one grouping variable.' );
end

if ngroup > 1 && ~anova
    warning ( 'Forcing the test to ANOVA to allow more than one grouping variable.' );
    anova = 1;
end

% If more than one grouping variable forbids permutations.
if ngroup > 1 && nperms > 0
    error ( 'Permutation test can not be performed with more than one grouping variable.' );
end


% Calculates the original p-value.
pval = stats ( data, group, covariates, anova, parametric );


% Checks if permutation statistics are requested.
if nperms
    
    % Initializes the permitation results.
    pperms = zeros ( size ( pval ) );
    
    % Iterates along permutations.
    for pindex = 1: nperms
        
        % Ramdomly re-labels the subjects.
        pgroup = group ( randperm ( numel ( group ) ) );
        
        % Calculates the permutation p-value.
        pperm  = stats ( data, pgroup, covariates, anova, parametric );
        
        % Stores the accumulated p-value.
        pperms = pperms + double ( pperm <= pval );
    end
    
    % Calculates the permutation-corrected p-value.
    pval = pperms / nperms;
end


% Function to get the statistics.
function pval = stats ( data, group, covariates, anova, parametric )

% Gets the number of grouping variables and of covariates, if any.
ngroup = size ( group, 2 );
ncov   = size ( covariates, 2 );

% Perfors ANOVA (with covariates) if required.
if anova
    
    % Initializes the input arguments.
    inputs    = { 'display', 'off' };
    
    % Rewrites the grouping variables to match anovan inputs.
    fullgroup = num2cell ( group, 1 );
    
    % If covariates, includes them as continuous grouping variables.
    if ncov
        fullgroup = cat ( 2, fullgroup, num2cell ( covariates, 1 ) );
        inputs    = cat ( 2, inputs, { 'continuous', ( 1: ncov ) + ngroup } );
    end
    
    % Gets the N-way ANOVA.
    [ pval, ~, stats ] = my_anovan ( data, fullgroup, inputs {:} );
    
    % Calculates the pairwise p-value with Tukey's HSD, if needed.
    pairwise  = my_multcompare ( stats, 'display', 'off' );
    if size ( pairwise, 1 ) == 1
        pval      = cat ( 1, pval, squeeze ( pairwise ( :, 6, : ) )' );
    else
        pval      = cat ( 1, pval, squeeze ( pairwise ( :, 6, : ) ) );
    end
    
% Otherwise performs a calculation for each pair of groups.
else
    
    % Gets the groups.
    groups  = unique ( group );
    ngroups = numel  ( groups );
    
    % Gets the shape of the comparison data.
    cshape  = size ( data );
    cshape (1) = 1;
    
    % Initializes the p-value variable.
    pval = [];
    
    % Goes through each pair or groups.
    for giindex = 1: ngroups
        for gjindex = giindex + 1: ngroups
            
            % Separates the groups.
            group1 = data ( group == giindex, : );
            group2 = data ( group == gjindex, : );
            
            % If parametric uses independent samples t-test.
            if parametric
                [ ~, cpval ] = ttest2 ( group1, group2 );
                
            % Otherwise uses a Wilcoxon's sum or ranks for each comparison.
            else
                cpval = myranksum ( group1, group2 );
            end
            
            % Adds the pair comparison to the matrix of p-values.
            pval = cat ( 1, pval, reshape ( cpval, cshape ) );
        end
    end
end
