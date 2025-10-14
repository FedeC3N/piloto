%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plot power spectrum of occipital sensors
%
% 01/10/2025
% Federico Ramirez-Torano
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
close all

% Adds the functions folders to PATH.
addpath(fullfile('..','SharedFunctions','functions'));
addpath(fullfile('..','SharedFunctions','mne_silent'));
addpath(fullfile('..','SharedFunctions','functions_eep'));
addpath(fullfile('..','SharedFunctions','fieldtrip-20250928'));
ft_defaults
clc

% Paths
config.path.raw         = fullfile('..','..','data','eeg','raw');
config.path.metadata    = fullfile('..','..','data','eeg','metadata');

% Action when the task have already been processed.
config.overwrite = false;

% Sets the list of component types.
% - Type 0 is always clean component ('Clean component').
% - Type 1 is always EOG component ('EOG component').
% - Type 2 is always EKG component ('EKG component').
config.comptypes = { 'Clean component' 'EOG component' 'EKG component' 'Noisy component' };

% Sets the visualization configuration parameters.
config.trialfun     = 'restingSegmentation';
config.segment      = 4;
config.overlap      = 0;
config.equal        = false;
config.padding      = 2;
config.addpadd      = true;

config.channel.groups = { 'EEG' };

% Sets the filter band.
config.filter.band    = [ 2 45 ];

% Sets the events to show (Inf for all).
config.event = Inf;

% Load the metadata
dataset_file = fullfile('..','..','data','eeg','dataset.mat');
load ( dataset_file );

% Filter the resting-state recordings
dummy           = regexp({dataset.task},'\w*EC');
resting_index   = ~cellfun(@isempty,dummy);
dataset         = dataset(resting_index);

% Goes through each subject and task.
n_freqs         = 176;
pow_spectrm_all = nan(n_freqs,numel(dataset));
for ifile = 1: numel ( dataset )
    
    % Get the current subject info
    current_folder          = dataset(ifile).folder;
    current_file            = dataset(ifile).file;
    current_subject         = dataset(ifile).subject;
    current_session         = dataset(ifile).session;
    current_task            = dataset(ifile).task;
    dummy_path_current_file = fullfile(current_folder,current_file);
    
    % Load metadata
    filename        = sprintf ('%s_%s_%s_metadata.mat',  ...
        current_subject, current_session, current_task );
    metadata_file   = fullfile(config.path.metadata,filename);
    metadata        = load (metadata_file);
    
    % Relevant info
    fileinfo    = metadata.fileinfo;
    header      = fileinfo.header;
    artinfo     = metadata.artinfo;
    
    % Load the recording
    cfg                 = [];
    cfg.dataset          = dummy_path_current_file;
    cfg.header           = header;
    cfg.channel          = {'eeg'};
    cfg.precision        = 'single';
    cfg.feedback         = 'no';
    wholedata            = ft_preprocessing ( cfg );
    
    % Segment in trials
    trialfun             = str2func ( config.trialfun );
    fileconfig           = config;
    fileconfig.dataset   = dummy_path_current_file;
    fileconfig.header    = wholedata.hdr;
    fileconfig.begtime   = fileinfo.begtime;
    fileconfig.endtime   = fileinfo.endtime;
    fileconfig.artifact  = artinfo.artifact;
    fileconfig.feedback  = 'no';
    fileconfig.trl       = trialfun ( fileconfig );
    
    trialdata            = ft_redefinetrial ( fileconfig, wholedata );
    clear wholedata;
    
    % Get the data matrix
    trialdata_matrix = cat(3,trialdata.trial{:});
    
    % Filter
    [b,a] = butter(4, config.filter.band/(header.Fs/2), 'bandpass');
    
    for itrial = 1 : size(trialdata_matrix,3)
        dummy                           =  filtfilt ( b,a, trialdata_matrix(:,:,itrial)' );
        trialdata_matrix(:,:,itrial)    = dummy';
    end
    
    % Remove padding
    padding_samples                                 = config.padding * header.Fs;
    trialdata_matrix(:,1:padding_samples-1,:)       = [];
    trialdata_matrix(:,end-padding_samples:end,:)   = [];
    
    % Reject badchannels and trials
    badchannel_index    = ismember(trialdata.label,metadata.chaninfo.bad);
    badtrial_index      = metadata.groupinfo.EEG.cleaninfo.trial.type ~= 0 ;
    trialdata_matrix    = trialdata_matrix(~badchannel_index,:,~badtrial_index);
    labels = trialdata.label(~badchannel_index);
    
    % Reject components
    mixing          = metadata.compinfo.SOBI.EEG.mixing;
    unmixing        = metadata.compinfo.SOBI.EEG.unmixing;
    comp_to_reject  = metadata.groupinfo.EEG.cleaninfo.comp.type ~= 0;
    
    dummy = reshape(trialdata_matrix,size(trialdata_matrix,1),[]);
    components = unmixing*dummy;
    components(comp_to_reject,:) = 0;
    dummy = mixing*components;
    cleandata_matrix = reshape(dummy,size(trialdata_matrix));
    
    % Estimate pow
    [num_channels,num_samples,num_trials]   = size(cleandata_matrix);
    win                                     = hann(num_samples,'periodic');         % column T x 1
    nfft                                    = 2^nextpow2(num_samples);
    num_freqs                               = floor(nfft/2)+1;
    Pxx_all                                 = zeros(num_channels, num_freqs, num_trials);
    f                                       = [];  % will be filled once
    for itrial = 1:num_trials
        for ichannel = 1:num_channels
            x = double(squeeze(cleandata_matrix(ichannel,:,itrial)));  % channel x trial signal
            [Pxx,f] = periodogram(x, win, nfft, header.Fs);
            Pxx_all(ichannel,:,itrial) = Pxx;           % store PSD
        end
    end
    f_index = (f >= 2) & (f <= 45);
    f    = f(f_index);                    
    Pxx_all  = Pxx_all(:, f_index, :); 
    
    % Average across trials and normalize across channels
    pow_spectrum = mean(Pxx_all,3);
    total_pow = sum(pow_spectrum,2);
    total_pow = repmat(total_pow,1,size(pow_spectrum,2));
    norm_pow_spectrum = pow_spectrum ./ total_pow;
    
    % Find the occipital channels
    dummy = regexp(labels,'\w*O\w*');
    occipital_index = ~cellfun(@isempty,dummy);
    
    % Save only the average pow_spectrum
    current_pow = nanmean(norm_pow_spectrum(occipital_index,:),1);
    pow_spectrm_all(:,ifile) = current_pow;
    
    
end

% Plot
figure
plot(f,pow_spectrm_all,'-')

for ifile = 1 : size(pow_spectrm_all,2)
    
    figure
    plot(f,pow_spectrm_all(:,ifile),'-')
    
end

