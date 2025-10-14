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
config.path.raw  = fullfile('..','..','data','eeg','raw');
config.path.metadata = fullfile('..','..','data','eeg','metadata');

% Sets the visualization configuration parameters.
config.trialfun       = 'restingSegmentation';
config.segment        = 20;
config.overlap        = 10;
config.equal          = false;
config.padding        = 2;
config.addpadd        = false;

config.channel.data   = { 'EEG' };
config.channel.ignore = {};

% Determines if the EEG data should be re-referenced.
config.channel.EEGref = 'average';
config.channel.hide   = 'zeros';

% Sets the filter band.
config.filter.band    = [ 2 45 ]; % [ 2 95 ];

% Action when the task have already been processed.
config.overwrite      = false;

% Load the metadata
dataset_file = fullfile('..','..','data','eeg','dataset.mat');
load ( dataset_file );

% Goes through each file.
for findex = 1: numel ( dataset )

    % Read the metadata
    current_folder = dataset(findex).folder;
    current_file = dataset(findex).file;
    current_subject = dataset(findex).subject;
    current_session = dataset(findex).session;
    current_task = dataset(findex).task;
    fprintf('Current file sub-%s_ses-%s_task-%s\n', current_subject,current_session,current_task);
    
    % Check if exist metadata file (mandatory)
    filename = sprintf ('%s_%s_%s_metadata.mat',  ...
        current_subject, current_session, current_task );
    metadata_file   = fullfile(config.path.metadata,filename);
    if ~exist(metadata_file)
        fprintf('Metadata file not created. Skip \n\n');
    end
    
    % Load metadata
    metadata = load (metadata_file);
    dummy_path_current_file = sprintf('%s/%s',current_folder,current_file);
    
    % Check if overwrite
    already_processed = strcmp(metadata.history.previous_steps,'Select badchannels');
    already_processed = sum(already_processed) > 0;
    if already_processed && ~config.overwrite
        fprintf('  Already processed. Do not overwrite\n\n');
        continue
    end
    
    fprintf ( 1, '  Reading data from disk.\n' );
    
    % Read the header
    header = ft_read_header(dummy_path_current_file);
    
    % Gets the data.
    cfg         = [];
    cfg.dataset = dummy_path_current_file;
    cfg.header  = header;
    
    wholedata   = my_read_data ( cfg );
    
    
    % Selects the channels.
    cfg                  = [];
    cfg.channel          = { 'EEG' 'EOG' 'ECG' };
    cfg.precision        = 'single';
    cfg.feedback         = 'no';
    
    wholedata            = ft_preprocessing ( cfg, wholedata );
    
    fprintf ( 1, '  Filtering the data in the band %0.0f - %0.0f Hz.\n', config.filter.band );
    
    % Calculates the optimal filter order from the desired padding.
    filtorder            = floor ( header.Fs * config.padding );
    downrate             = floor ( header.Fs / ( 2 * config.filter.band ( end ) ) );
    
    % Filters and downsamples the data.
    fir                  = fir1 ( filtorder, config.filter.band / ( wholedata.fsample / 2 ) );
    wholedata            = myft_filtfilt ( fir, 1, wholedata );
    wholedata            = my_downsample ( wholedata, downrate );
     
    % Extracts the overlapping epochs for the artifact revision.
    trialfun             = str2func ( config.trialfun );
    
    fileconfig           = config;
    fileconfig.dataset   = dummy_path_current_file;
    fileconfig.header    = wholedata.hdr;
    fileconfig.begtime   = NaN;
    fileconfig.endtime   = NaN;
    fileconfig.feedback  = 'no';
    
    fileconfig.trl       = trialfun ( fileconfig );
    
    trialdata            = ft_redefinetrial ( fileconfig, wholedata );
    clear wholedata;
    
    fprintf ( 1, '  Displaying the data.\n' );
    
    % Gets the list of non-MEG channels (EOG, EKG, etc.).
    MEG                 = find ( ismember ( trialdata.label', ft_channelselection ( { 'MEG' },             trialdata.label ) ) );
    MEGMAG              = find ( ismember ( trialdata.label', ft_channelselection ( { 'MEGMAG' },          trialdata.label ) ) );
    MEGGRAD             = find ( ismember ( trialdata.label', ft_channelselection ( { 'MEGGRAD' },         trialdata.label ) ) );
    EEG                 = find ( ismember ( trialdata.label', ft_channelselection ( { 'EEG' },             trialdata.label ) ) );
    physio              = find ( ismember ( trialdata.label', ft_channelselection ( { 'EOG' 'ECG' 'EMG' }, trialdata.label ) ) );
    
    % Sets the colors according to the channel type.
    chancol             = ones ( numel ( trialdata.label ), 1 );
    chancol ( MEG )     = 1;
    chancol ( MEGMAG )  = 2;
    chancol ( MEGGRAD ) = 3;
    chancol ( EEG )     = 4;
    chancol ( physio )  = 5;
    
    % Displays the EEG data to mark badchannels
    cfg = [];
    cfg.channel         = 1: 30;
    cfg.physio          = physio;
    cfg.eegscale        = 1e-13;
    cfg.eogscale        = 1e-13;
    cfg.ecgscale        = 1e-13;
    cfg.ylim            = [ -1 1 ] * 4e-12;
    cfg.plotlabels      = 'yes';
    cfg.viewmode        = 'vertical';
    cfg.continous       = 'yes';
    cfg.colorgroups     = chancol;
    cfg.channelcolormap = [   0 114 189;   0 114 189; 162  20  47; 217  83  25; 126  47 142 ] / 255;
    cfg.ploteventlabels = 'colorvalue';
    cfg.selectmode      = 'markbadchannel';
    
    cfg                 = my_databrowser ( cfg, trialdata );
    
    % Save the metadata
    % artinfo
    metadata.chaninfo.bad = cfg.badchan;
    
    % history
    metadata.history.last_step = 'Select badchannels';
    metadata.history.previous_steps(end+1) = {metadata.history.last_step};
    metadata.history.date = datestr ( now );
    
    % Saves the output data.
    save ( '-v6', metadata_file, '-struct', 'metadata' );
    fprintf('\n\n')
    
end
