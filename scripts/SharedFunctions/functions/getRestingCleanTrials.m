function trl = getRestingCleanTrials ( metadata, artifacts, header )

%getRestingCleanTrials Gets resting trials avoiding artifacts.
%   
%   Usage: getRestingCleanTrials (METADATA,ARTIFACTS)
%   
%   This function obtains a list of not necesarilly consecutive and
%   non-artifacted trials of duration t_frag between the times t_on and
%   t_off from the file fif_path.
%   
%   t_frag, t_on, t_off and fif_path are fields of the structure METADATA.
%   ARTIFACTS is a n x 2 matrix with the beggining and ending sample of
%   the artifact detected between t_on and t_off.

% Initializes the trl matrix.
trl       = zeros ( 0, 2 );

% Gets the header to obtain the sampling rate.
if nargin < 3, header    = ft_read_header ( metadata.dataset ); end
if nargin < 2, artifacts = zeros ( 0, 2 );                      end

% Gets the samples limit.
lims      = round ( [ metadata.start metadata.end ] * header.Fs );

% Checks that the limits are valid.
if lims (1) < 1
    warning ( 'Given start time is before the beginning of the data.' );
    lims (1) = 1;
end
if lims (2) > header.nSamples
    warning ( 'Given end time is after the ending of the data.' );
    lims (2) = header.nSamples;
end

% Creates a vector with the resting samples.
span      = lims (1): lims (2);

% Concatenates all the artifact matrices.
if isstruct ( artifacts )
    artifacts = struct2cell ( artifacts );
end
if isstruct ( artifacts {1} )
    artifacts = cellfun ( @(artifacts) artifacts.artifact, artifacts, 'UniformOutput', false );
end
if iscell   ( artifacts )
    artifacts = cat ( 1, artifacts {:} );
end

% Replaces the samples for NaNs when an artifact is present.
for artifact = 1: size ( artifacts, 1 )
    span ( span >= artifacts ( artifact, 1 ) & span <= artifacts ( artifact, 2 ) ) = NaN;
end

% Gets the list of artifact endings.
transitions = find ( cat ( 2, 1, diff ( ~isnan ( span ) ) ) > 0 );

% Goes through all the artifact ends.
for first = transitions
    
    % Infinite loop :D
    while true
        
        % Sets the local span.
        last  = first + metadata.segment * header.Fs - 1;
        
        % If the segment contains an artifact exits.
        if last > numel ( span ) || any ( isnan ( span ( first: last ) ) ), break, end
        
        % Stores the trial.
        trl   = cat ( 1, trl, [ first last ] );
        
        % Moves the pointer.
        first = last + 1; %#ok<FXSET>
    end
end

% Replaces the timespan indexes to samples and adds the baseline duration.
trl = span ( trl );
trl ( end, 3 ) = 0;