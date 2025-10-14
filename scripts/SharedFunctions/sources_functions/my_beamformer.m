function sources = my_beamformer ( cfg, grid, data )

if ~ft_datatype ( data, 'timelock' )
    error ( 'This function only works with LCMV beam former for now.' )
end

% Generates boolean flags.
project   = isfield ( cfg, 'projectmom' )   && strcmp ( cfg.projectmom,   'yes' );
keepnoise = isfield ( cfg, 'projectnoise' ) && strcmp ( cfg.projectnoise, 'yes' );
keepcov   = isfield ( cfg, 'keepcov' )      && strcmp ( cfg.keepcov,      'yes' );
powsvd    = isfield ( cfg, 'powmethod' )    && strcmp ( cfg.powmethod,    'lambda1' );


% Gets the whitener, if provided.
if isfield ( cfg, 'subspace' ) && ~isempty ( cfg.subspace )
    whitener = cfg.subspace;
else
    whitener = eye ( size ( data.cov ) );
end

% Calculates the whitened covariance matrix.
wcov = whitener * data.cov * whitener';


% Gets the regularization parameter, if provided.
if isfield ( cfg, 'lambda' )
    
    % If numeric, takes the value of lambda.
    if isnumeric ( cfg.lambda )
        lambda = cfg.lambda;
        
    % If percent, calculates the value of lambda.
    elseif ischar ( cfg.lambda ) && strcmp ( cfg.lambda ( end ), '%' )
        ratio  = str2double ( cfg.lambda ( 1: end - 1 ) ) / 100;
        lambda = ratio * trace ( wcov ) / size ( wcov, 1 );
    end
else
    lambda = 0;
end

if ~isequal ( whitener, diag ( diag ( whitener ) ) ) && lambda ~= 0
    warning ( 'Using both data whitening and Tikhonov regularization.' )
end

% Calculates the (regularized) inverse of the covariance matrix.
% icov = pinv ( whitener * data.cov * whitener' );
icov = pinv ( wcov + lambda * eye ( size ( wcov ) ) );


% Estimates the channel-level noise.
noise = svd ( data.cov );
noise = noise ( end );
noise = max ( noise, lambda );


% Initializes the sources structure.
sources          = [];
sources.label    = grid.label;
sources.pos      = grid.pos;
sources.inside   = grid.inside;
sources.filter   = cell ( size ( grid.leadfield (:) ) );
sources.pow      = nan  ( size ( grid.leadfield (:) ) );

if keepnoise
    sources.noise    = nan  ( size ( grid.leadfield (:) ) );
end
if keepcov
    sources.cov      = cell ( size ( grid.leadfield (:) ) );
end
if keepnoise && keepcov
    sources.noisecov = cell ( size ( grid.leadfield (:) ) );
end

% Goes through each source.
for c = find ( grid.inside (:) )'
    
    % Gets the (whitened) leadfield for the current dipole.
    dipleadfield = whitener * grid.leadfield { c };
    
    % Calculates the (unwhitened) beam former filter.
    dipfilter    = pinv ( dipleadfield' * icov * dipleadfield ) * dipleadfield' * icov * whitener;
    
    % Calculates the sources' power and noise.
    dippowercov  = dipfilter * data.cov * dipfilter';
    dipnoisecov  = noise * ( dipfilter * dipfilter' );
    
    % Projects over the dominant direction, if requested.
    if project
        [ u, ~ ]     = svd ( dippowercov );
        dipfilter    = u ( :, 1 )' * dipfilter;
        dippowercov  = u ( :, 1 )' * dippowercov * u ( :, 1 );
        dipnoisecov  = u ( :, 1 )' * dipnoisecov * u ( :, 1 );
    end
    
    % Calculates the dipole power.
    if powsvd
        dippow       = max ( svd ( dippowercov ) );
        dipnoise     = max ( svd ( dipnoisecov ) );
    else
        dippow       = trace ( dippowercov );
        dipnoise     = trace ( dipnoisecov );
    end
    
    % Stores the mandatory source-space information.
    sources.filter   { c } = dipfilter;
    sources.pow      ( c ) = dippow;
    
    % Stores the extra information.
    if keepnoise
        sources.noise    ( c ) = dipnoise;
    end
    if keepcov
        sources.cov      { c } = dippowercov;
    end
    if keepnoise && keepcov
        sources.noisecov { c } = dipnoisecov;
    end
end
