function trl = trialfunResting ( cfg )
% trialfunResting ( config )
%
% Segments the data according to the 'config' structure:
% - config.start:     Beggining time.
% - config.end:       Ending time.
% - config.segment:   Padding time.
% - config.trllen:    Trial duration.
% - config.overlap:   Trials overlap.
% - config.eqtrl:     If true, all the trials must have the same length.
%
% This function segments the data in continous (non-)overlapping segments
% of 'segment' length, the first segment starting in 'start' adn ending in
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

% Gets the file header.
hdr   = ft_read_header ( cfg.dataset );

% Stablish the time of interest.
tin   = cfg.start;
tf    = cfg.end;
tpad  = cfg.padding;
tlen  = cfg.segment;
tstep = cfg.segment - cfg.overlap;
eqtrl = cfg.eqtrl;

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

% Fills the trials.
for i = 1: ntrl
    begsample    = sin + ( i - 1 ) * sstep + 1;
    endsample    = sin + ( i - 1 ) * sstep + slen;
    begsample    = begsample - spad;
    endsample    = endsample + spad;
    trl ( i, : ) = [ begsample endsample off ];
end

% Removes the trials beyond the end of the data.
ltrl  = find ( trl ( :, 2 ) > sf + spad, 1, 'first' );
trl   ( ltrl + 1: end, : ) = [];

% Makes sure the last trial doesn't go beyond the data.
if eqtrl
    trl ( end, 1  ) = sf - slen - spad;
    trl ( end, 2  ) = sf + spad;
else
    trl ( end, 2  ) = min ( trl ( end, 2  ), sf + spad );
end
