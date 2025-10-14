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

% Action when the task have already been processed.
config.overwrite      = true;

% Sets the segmentation options.
config.trialfun       = 'restingSegmentation';
config.segment        = 4.0;
config.padding        = 2.0;
config.overlap        = 0.0;

% Sets the sketch building options.
config.freqband       = [ 2 45 ];

config.channel.groups = { 'EEG' };
config.channel.ignore = {};

% Load the metadata
dataset_file = fullfile('..','..','data','eeg','dataset.mat');
load ( dataset_file );
% Goes through each subject and task.
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
    dummy_path_current_file = fullfile(current_folder,current_file);
    
    % Check if overwrite
    already_processed = strcmp(metadata.history.previous_steps,'Sketch Extraction');
    already_processed = sum(already_processed) > 0;
    if already_processed && ~config.overwrite
        fprintf('  Already processed. Do not overwrite\n\n');
        continue
    end
    
    % Checks that the selected channel group is present in the SOBI data.
    if ~any ( ismember ( config.channel.groups, fieldnames ( metadata.compinfo.SOBI ) ) )
        fprintf ( 1, 'Ignoring %s (no SOBI information for the selected channel groups).\n', msgtext );
        continue
    end
    
    % Gets the files length as a vector.
    headers             = [ metadata.fileinfo.header ];
    samples             = [ headers.nSamples ];
    
    % Reserves memory for the trialdata.
    trialdefs           = cell ( numel ( metadata.fileinfo ), 1 );
    trialfiles          = cell ( numel ( metadata.fileinfo ), 1 );
    trialpads           = cell ( numel ( metadata.fileinfo ), 1 );
    trialtypes          = cell ( numel ( metadata.fileinfo ), 1 );
    erfdatas            = cell ( numel ( metadata.fileinfo ), 1 );
    freqdatas           = cell ( numel ( metadata.fileinfo ), 1 );
        
    % Loads the current epochs and artifact definitions.
    fileinfo             = metadata.fileinfo ;
    artinfo              = metadata.artinfo;
    
    fprintf ( 1, '    Reading data from disk.\n' );
    
    % Gets the MEG data.
    cfg                   = [];
    cfg.dataset           = dummy_path_current_file;
    cfg.header            = fileinfo.header;
    
    wholedata             = my_read_data ( cfg );
    
    % Selects the channels.
    cfg                  = [];
    badchannels          = strcat('-',metadata.chaninfo.bad);
    desired_channel      = cat(2,{'eeg'},badchannels);
    cfg.channel          = desired_channel;
    cfg.precision        = 'single';
    cfg.feedback         = 'no';
    
    wholedata            = ft_preprocessing ( cfg, wholedata );
    
    
    fprintf ( 1, '    Generating the trial definition. ' );
    
    % Extracts the trials.
    trialfun              = str2func ( config.trialfun );
    
    fileconfig            = config;
    fileconfig.dataset    = dummy_path_current_file;
    fileconfig.header     = fileinfo.header;
    fileconfig.event      = fileinfo.event;
    fileconfig.begtime    = fileinfo.begtime;
    fileconfig.endtime    = fileinfo.endtime;
    fileconfig.artifact   = artinfo.artifact;
    fileconfig.feedback   = 'no';
    
    trialdef              = trialfun ( fileconfig );
    
    % If no trials ignores the file.
    if isempty ( trialdef )
        fprintf ( 1, 'No trials detected. Skipping file.\n' );
        continue
    end
    
    % Initializes the file configuration structure.
    fileconfig            = config;
    
    fprintf ( 1, '%i trials found.\n', size ( trialdef, 1 ) );
    
    
    % Gets the defined trials.
    fileconfig.feedback   = 'no';
    fileconfig.trl        = trialdef;
    
    freqdata              = ft_redefinetrial ( fileconfig, wholedata );
    clear wholedata
    
    
    % Creates a dummy structure with no data.
    erfdata               = freqdata;
    erfdata.trial         = cellfun ( @(x) zeros ( cat ( 2, size ( x ), 0 ) ), erfdata.trial, 'UniformOutput', false );
    erfdata.trialdimord   = '{rpt}_chan_time';
    
    % Removes the 'cfg' field.
    erfdata               = rmfield ( erfdata,  'cfg' );
    
    
    % Calculates the spectra of the data.
    cfg                   = [];
    cfg.method            = 'mtmfft';
    cfg.taper             = 'hamming';
    cfg.output            = 'fourier';
    cfg.foilim            = config.freqband;
    cfg.keeptrials        = 'yes';
    cfg.feedback          = 'no';
    
    freqdata              = ft_freqanalysis ( cfg, freqdata );
    
    % Removes the 'cfg' field.
    freqdata              = rmfield ( freqdata,  'cfg' );
    
    % Rewrites the data as 'single' to save space.
    freqdata.fourierspctrm = single ( freqdata.fourierspctrm );
    
    
    % Writes the trial's metadata.
    trialtype             = zeros ( numel ( erfdata.trial ), 1 );
    trialpad              = zeros ( numel ( erfdata.trial ), 1 );
    trialfile             = repmat ( findex,   numel ( erfdata.trial ), 1 );

    % Generates the trial information structure from the trial data.
    trialinfo             = [];
    trialinfo.trialdef    = trialdef;
    trialinfo.trialfile   = trialfile;
    trialinfo.trialpad    = trialpad;
    
    % Initilizes the clearn trials and components structure.
    cleaninfo             = [];
    cleaninfo.trial.types = { 'Clean trial' 'Noisy trial' };
    cleaninfo.trial.type  = trialtype;
    
    
    % Goes through each channel group.
    for chindex = 1: numel ( config.channel.groups )
        
        channel             = config.channel.groups { chindex };
        
        if ~ismember ( channel, fieldnames ( metadata.compinfo.SOBI ) )
            fprintf ( 1, '  Ignoring channel group ''%s'' (no SOBI information).\n', channel );
            continue
        end
        
        fprintf ( 1, '  Saving channel group ''%s''.\n', channel );
        
        % Keeps only the selected channel group data.
        cfg                   = [];
        cfg.channel           = cat ( 2, { channel }, { 'EOG' 'ECG' 'EMG' } );
        cfg.feedback          = 'no';
        
        grouperf              = ft_selectdata ( cfg, erfdata   );
        groupfreq             = ft_selectdata ( cfg, freqdata  );
        
        % Removes the 'cfg' field.
        grouperf              = rmfield ( grouperf,   'cfg' );
        groupfreq             = rmfield ( groupfreq,  'cfg' );
        
        % Keeps only the selected channel group SOBI information.
        compinfo              = metadata.compinfo;
        compinfo.SOBI         = metadata.compinfo.SOBI.( channel );
        
        % Adds the clean components to the information structure.
        cleaninfo.comp.types  = compinfo.types;
        cleaninfo.comp.type   = compinfo.SOBI.type;
        
        % Fills the group information.
        groupinfo             = [];
        groupinfo.(channel).erfdata     = grouperf;
        groupinfo.(channel).freqdata    = groupfreq;
        groupinfo.(channel).trialinfo   = trialinfo;
        groupinfo.(channel).cleaninfo   = cleaninfo;
        
    end
    
    % Add the info
    metadata.groupinfo = groupinfo;

    % Update the history
    metadata.history.last_step = 'Sketch Extraction';
    metadata.history.previous_steps(end+1) = {metadata.history.last_step};
    metadata.history.date = datestr ( now );
    
    % Saves the current group epoch data.
    save ( '-v6', metadata_file, '-struct', 'metadata' );
    
end
