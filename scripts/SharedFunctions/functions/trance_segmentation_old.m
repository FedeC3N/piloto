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
% 1 baseline of 10 seconds
% 6 trials of 10 seconds
epochs   = zeros ( 7 , 3 );

% Get the 10 seconds epochs
epochs_beggining = (0:6)*epoch_length + 1;
epochs_end = epochs_beggining + epoch_length - 1;

% Get the 9 second epoch
epochs(:,1) = epochs_beggining;
epochs(:,2) = epochs_end;

