function trl = trialfunResting ( cfg, hdr )
% trialfunResting ( config )
%
% Segments the data according to the 'config' structure:
% - config.dataset:   MEG dataset.
% - config.start:     Beggining time.
% - config.end:       Ending time.
% - config.segment:   Padding time.
% - config.trllen:    Trial duration.
% - config.overlap:   Trials overlap.
% - config.eqtrl:     If true, all the trials must have the same length.
% - config.artifact:  Artifacts definition.
%
% This function segments the data in continous (non-)overlapping segments
% of 'segment' length, the first segment starting in 'start' and ending in
% 'end'.
% If 'padding' is non-zero, segments are expanded by 'padding' seconds, and
% the offset field in the 'trl' matrix is set accordingly.
%
% - If 'start' is smaller than 'padding' its value is set to 'padding' + 1.
% - If 'end' is greater than data length - 'padding', its value is set to
%   data length - 'padding' - 1.
%
% 'eqtrl' defines the behaviour to handle the last incomplete trial.
% - If 'eqtrl' is 'true', the begining of the last trial is modified to
%   have a complete trial ending in 'end' + 'padding'.
% - If 'eqtrl' is 'false', the last trial is incomplete and ends in 'end' +
%   'padding'.

% Fields 'dataset' and 'segment' are mandatory.
if ~isfield ( cfg, 'dataset'  ), error ( '''dataset'' field is mandatory.' ); end
if ~isfield ( cfg, 'segment'  ), error ( '''segment'' field is mandatory.' ); end

% Fulfills the configuration structure.
if ~isfield ( cfg, 'start'    ), cfg.start    = 0;    end
if ~isfield ( cfg, 'end'      ), cfg.end      = inf;  end
if ~isfield ( cfg, 'outpad'   ), cfg.outpad   = 0;    end
if ~isfield ( cfg, 'padding'  ), cfg.padding  = 0;    end
if ~isfield ( cfg, 'overlap'  ), cfg.overlap  = 0;    end
if ~isfield ( cfg, 'eqtrl'    ), cfg.eqtrl    = true; end
if ~isfield ( cfg, 'artifact' ), cfg.artifact = [];   end

% Gets the file header, if not supplied.
if nargin < 2, hdr = ft_read_header ( cfg.dataset ); end

% Gets the configuration.
tin   = cfg.start;
tf    = cfg.end;
tpad  = cfg.padding + cfg.outpad;
tlen  = cfg.segment;
tstep = cfg.segment - cfg.overlap;
eqtrl = cfg.eqtrl;
sarts = cfg.artifact;

% Gets the samples related to those times.
sin   = round ( tin   * hdr.Fs );
sf    = round ( tf    * hdr.Fs );
spad  = round ( tpad  * hdr.Fs );
slen  = round ( tlen  * hdr.Fs );
sstep = round ( tstep * hdr.Fs );

if ( sin < spad )
    warning ( 'myApp:argChk', 'The initial time is smaller than the padding.' )
    sin = spad;
end

if ( sf  > hdr.nSamples - spad - 1 )
    warning ( 'myApp:argChk', 'The final time is greater than the data.' )
    sf  = hdr.nSamples - spad - 1;
end

% Gets the maximal number of (overlapping) trials.
ntrl  = ceil ( ( ( sf + spad ) - sin - ( slen - sstep ) ) / sstep );
off   = -spad;

% Reserves memory for the trial information.
trl   = zeros ( ntrl, 3 );


% Creates the artifacts vector.
span = 1: sf + spad;
span ( sin ) = NaN;

% Fills the artifacts vector, if the artifact definition was provided.
if numel ( sarts )
    
    % Concatenates all the artifact matrices.
    if isstruct ( sarts )
        sarts = struct2cell ( sarts );
    end
    if isstruct ( sarts {1} )
        sarts = cellfun ( @(artifacts) artifacts.artifact, sarts, 'UniformOutput', false );
    end
    if iscell   ( sarts )
        sarts = cat ( 1, sarts {:} );
    end
    
    % Replaces the samples for NaNs when an artifact is present.
    for artifact = 1: size ( sarts, 1 )
        span ( sarts ( artifact, 1 ): sarts ( artifact, 2 ) ) = NaN;
    end
    
    % Marks as artifacts the data outside sin-sf.
    span ( 1:  sin ) = NaN;
    span ( sf: end ) = NaN;
end

% Gets the list of artifact endings.
begs  = find ( diff ( ~isnan ( span ) ) > 0 );

% Removes the gaps that can't fit a segment.
begs  ( diff ( begs ) < slen ) = [];


% Initializes the segment index.
segment = 1;

% Goes through the artifact ends.
for beg = begs
    
    for i = 1: ntrl
        
        % Sets the segment edges.
        begsample = beg + ( i - 1 ) * sstep + 1;
        endsample = beg + ( i - 1 ) * sstep + slen;
        
        % If the segment contains artifacts, jumps to the next gap.
        maxlen    = find ( isnan ( span ( begsample: end ) ), 1, 'first' );
        if maxlen < slen, break, end
        
        % Adds the padding and saves the segment.
        begsample = begsample - spad;
        endsample = endsample + spad;
        
        trl ( segment, : ) = [ begsample endsample off ];
        
        % Updates the segment index.
        segment   = segment + 1;
    end
end

% Removes the empty trials.
ltrl  = find ( any ( trl, 2 ), 1, 'last' );
trl   ( ltrl + 1: end, : ) = [];

% Removes the trials beyond the end of the data.
ltrl  = find ( trl ( :, 2 ) > sf + spad, 1, 'first' );
trl   ( ltrl + 1: end, : ) = [];

% Makes sure the last trial doesn't go beyond the data.
if trl ( end, 2 ) > sf + spad
    if eqtrl
        trl ( end, 1  ) = sf - slen - spad;
        trl ( end, 2  ) = sf + spad;
    else
        trl ( end, 2  ) = min ( trl ( end, 2  ), sf + spad );
    end
end

% Removes the outter padding, if any.
trl = bsxfun ( @minus, trl, [ -1 1 -1 ] * cfg.outpad * hdr.Fs );
