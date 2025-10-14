function epochs = trance_segmentation ( cfg )
% Segments the data according to the 'config' structure.
% Structure fields are:
% - cfg.dataset:  Data file.
% - cfg.epoch:  Epoch duration, in seconds.
% - cfg.padding:  Padding in seconds. It will be added equal on both sides.
% It will not be removed at the end of the trial definition.
%
% This function segments the data in consecutive
% (non)overlapping segments of 'segment' seconds length, the first segment
% trying to start in 'begtime'. The segments can not extend beyond
% 'endtime'. All the segments must be surrounded by, at least,
% 'padding' seconds of data.
%
% The output is a Nx3 matrix, indicating in each row the first sample, the
% last sample and the number of pre-zero samples of a given segment.
%

% 'dataset' and 'segment' are mandatory fields.
if ~isfield ( cfg, 'dataset' ), error ( 'No dataset provided,' );        end
if ~isfield ( cfg, 'epoch' ), error ( 'No segment length provided,' ); end

% Read the header
header = ft_read_header(cfg.dataset);
Fs = header.Fs;

% Gets the samples related to those times.
epoch_length = cfg.epoch * Fs;

% Reserves memory for the epochs.
% If the last trial is not big enough, discard it
num_epochs = floor(cfg.song_length/cfg.epoch);
epochs   = zeros ( num_epochs + 1 , 3 ); % baseline + num_epochs

% Get the epochs beggining and end
% Remember the baseline is 10.000 samples always (10 seconds)
epochs_beggining = 10000 + (0:num_epochs-1)*epoch_length + 1;
epochs_end = epochs_beggining + epoch_length - 1;

% Store it in a fieldtrip way
% Baseline
epochs(1,1) = 1;
epochs(1,2) = 10000;
% Song
epochs(2:end,1) = epochs_beggining;
epochs(2:end,2) = epochs_end;




