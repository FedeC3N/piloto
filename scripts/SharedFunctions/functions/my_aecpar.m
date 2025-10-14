function aec = my_aec ( data, sampok, ortho, nsmooth )
% my_aec Calculation of Amplitude Envelope Correlation (AEC)
%
% Use:
%  aec = my_aec ( data, cleansamples, orthogonalize, smooth_window )
%  aec = my_aec ( data, cleansamples, orthogonalize, smooth_window )
% 
% * data is a nsamples x nsignals matrix.
% * cleansamples is a nsamples x 1 binary mask for the clean samples.
% * ortogonalize is a binary flag indicating if the data should be
%   (pairwise) orthogonalized prior to calculation of AEC. Default is to
%   orthogonalize.
% * smooth_window is the window length used for the smoothing. If 0, no
%   smooth is performed. Default is 250 samples.

if nargin < 2, error ( 'Not enough input arguments.' ); end
if nargin < 3, ortho = true; end
if nargin < 4, nsmooth = 250; end

% Gets the size of the data.
nsamp = size ( data, 1 ); %#ok<NASGU>
nsign = size ( data, 2 );
nsok  = sum ( sampok ); %#ok<NASGU>

% Demeans the data.
data  = bsxfun ( @minus, data, mean ( data, 1 ) );

% Gets the Hilbert transform.
xs    = hilbert ( data );

% Gets the envelope.
xes   = abs ( xs );

% Smooths the envelope, if requested.
if nsmooth
    xes   = my_masmooth ( xes, nsmooth );
end

% Removes the invalid points.
xes   = xes ( sampok, : );


% Calculates the AEC without smoothing, if requested.
if ~ortho
    aec   = corr ( xes );
    return
end


% Calculates the projection (betas) for each pair of signals.
projs = data ( sampok, : )' * data ( sampok, : );
betas = bsxfun ( @rdivide, projs, diag ( projs ) );

% Reserves memory for the correlation matrix.
aec   = zeros ( nsign );

% Goes through each regressed source.
parfor ireg = 1: nsign
    
    % Removes the current signal from all the signals.
    yrs   = xs - xs ( :, ireg ) * betas ( ireg, : );
    
    % Gets the envelope.
    yres  = abs ( yrs );
    
    % Smooths the envelope, if requested.
    if nsmooth
        yres  = my_masmooth ( yres, nsmooth );
    end
    
    % Removes the invalid points.
    yres  = yres ( sampok, : );
    
    % Calculates the correlation.
    aec ( ireg, : ) = corr ( xes ( :, ireg ), yres );
end

% Simetrices the correlation matrix.
aec   = ( aec + aec' ) ./2;
