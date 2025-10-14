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
config.overwrite      = false;

% Sets the visualization configuration parameters.
config.trialfun       = 'restingSegmentation';
config.segment        = 20;
config.overlap        = 10;
config.equal          = false;
config.padding        = 2;
config.addpadd        = false;

config.channel.data   = {'EEG' };
config.channel.ignore = {};

% Sets the filter band.
config.filter.band    = [ 2 45 ];

% Sets the events to show (Inf for all).
config.event          = Inf;

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
    already_processed = strcmp(metadata.history.previous_steps,'Artifact revision');
    already_processed = sum(already_processed) > 0;
    if already_processed && ~config.overwrite
        fprintf('  Already processed. Do not overwrite\n\n');
        continue
    end
    
    fprintf ( 1, 'Working on %s.\n', metadata.fileinfo.file );
    
    % Gets the files length as a vector.
    headers             = [ metadata.fileinfo.header ];
    samples             = [ headers.nSamples ];
    
    % Calculates the optimal downsampling rate for the frequency band.
    downrate            = floor ( headers (1).Fs / ( 2 * config.filter.band ( end ) ) );
    
    % Calculates the optimal filter order from the desired padding.
    filtorder           = floor ( headers (1).Fs * config.padding );
    
    % Gets the current epochs and artifact definitions.
    fileinfo             = metadata.fileinfo;
    artinfo              = metadata.artinfo;
    
    fprintf ( 1, '    Reading data from disk.\n' );
    
    % Gets the MEG data.
    cfg                  = [];
    cfg.dataset          = dummy_path_current_file;
    cfg.header           = fileinfo.header;
    
    wholedata            = my_read_data ( cfg );
    
    % Selects the channels.
    cfg                  = [];
    badchannels          = strcat('-',metadata.chaninfo.bad);
    desired_channel      = cat(2,{'eeg'},badchannels);
    cfg.channel          = desired_channel;
    cfg.precision        = 'single';
    cfg.feedback         = 'no';
    
    wholedata            = ft_preprocessing ( cfg, wholedata );
    
    
    fprintf ( 1, '    Filtering the data in the band %0.0f - %0.0f Hz.\n', config.filter.band );
    
    % Filters and downsamples the data.
    fir                  = fir1 ( filtorder, config.filter.band / ( wholedata.fsample / 2 ) );
    wholedata            = myft_filtfilt ( fir, 1, wholedata );
    wholedata            = my_downsample ( wholedata, downrate );
    
    
    % Extracts the overlapping epochs for the artifact revision.
    trialfun             = str2func ( config.trialfun );
    
    fileconfig           = config;
    fileconfig.dataset   = dummy_path_current_file;
    fileconfig.header    = wholedata.hdr;
    fileconfig.begtime   = fileinfo.begtime;
    fileconfig.endtime   = fileinfo.endtime;
    fileconfig.feedback  = 'no';
    
    fileconfig.channel.bad = metadata.chaninfo.bad;
    
    fileconfig.trl       = trialfun ( fileconfig );
    
    trialdata            = ft_redefinetrial ( fileconfig, wholedata );
    clear wholedata;
    
    % Removes the 'cfg' field.
    trialdata            = rmfield ( trialdata, 'cfg' );
    trialdata.trial      = cellfun ( @single, trialdata.trial, 'UniformOutput', false );
    
    % Resamples and corrects the events definitions.
    event                = fileinfo.event (:);
    
    if numel ( event )
        [ event.sample ]      = my_deal ( ceil ( [ event.sample ] / downrate ) );
    end
    
    % Resamples and corrects the artifact definitions.
    artifact              = artinfo.artifact;
    arttypes              = fieldnames ( artifact );
    
    for aindex = 1: numel ( arttypes )
        arttype                       = arttypes { aindex };
        artifact.( arttype ).artifact = ceil ( artifact.( arttype ).artifact / downrate );
    end
    
    % Removes the undesired events, if requested.
    if config.event < Inf
        event               = event ( ismember ( [ event.value ], config.event ) );
    end
      
    fprintf ( 1, '  Combining the components for all channel groups.\n' );
    
    % Selects the SOBI components for the last step.
    SOBI                = metadata.compinfo.SOBI;
    
    % Combines the components for all the chanel types.
    channels            = fieldnames ( SOBI );
    comptypes           = cell ( numel ( channels ), 1 );
    compposs            = cell ( numel ( channels ), 1 );
    topolabels          = cell ( numel ( channels ), 1 );
    unmixings           = cell ( numel ( channels ), 1 );
    mixings             = cell ( numel ( channels ), 1 );
    EOGcomps            = cell ( numel ( channels ), 1 );
    EKGcomps            = cell ( numel ( channels ), 1 );
    
    % Goes through each channel type.
    for chindex = 1: numel ( channels )
        
        channel                = channels { chindex };
        
        % Gets the values for the current channel type.
        comptypes  { chindex } = chindex * ones ( numel ( SOBI.( channel ).topolabel ), 1 );
        compposs   { chindex } = ( 1: numel ( SOBI.( channel ).topolabel ) )';
        topolabels { chindex } = SOBI.( channel ).topolabel;
        mixings    { chindex } = SOBI.( channel ).mixing;
        unmixings  { chindex } = SOBI.( channel ).unmixing;
        
        % Gets the values for the rejected components.
        if isfield ( SOBI.( channel ), 'type' )
            EOGcomps { chindex } = find ( SOBI.( channel ).type == 1 );
            EOGcomps { chindex } = cat ( 2, chindex * ones ( size ( EOGcomps { chindex } ) ), EOGcomps { chindex } );
            
            EKGcomps { chindex } = find ( SOBI.( channel ).type == 2 );
            EKGcomps { chindex } = cat ( 2, chindex * ones ( size ( EKGcomps { chindex } ) ), EKGcomps { chindex } );
        else
            EOGcomps { chindex } = zeros ( 0, 2 );
            EKGcomps { chindex } = zeros ( 0, 2 );
        end
    end
    
    % Concatenates the matrices.
    comptype            = cat     ( 1, comptypes  {:} );
    comppos             = cat     ( 1, compposs   {:} );
    topolabel           = cat     ( 1, topolabels {:} );
    mixing              = blkdiag ( mixings       {:} );
    unmixing            = blkdiag ( unmixings     {:} );
    EOGcomp             = cat     ( 1, EOGcomps   {:} );
    EKGcomp             = cat     ( 1, EKGcomps   {:} );
    
    % Tries to sort the components.
    %%% TO DO
    
    
    fprintf ( 1, '  Removing the artifact components.\n' );
    
    % Calculates the position of the current EOG and EKG components.
    EOGcomp             = find ( ismember ( [ comptype comppos ], EOGcomp, 'rows' ) );
    EKGcomp             = find ( ismember ( [ comptype comppos ], EKGcomp, 'rows' ) );
    
    % Gets the components for all the channel types at once.
    cfg = [];
    cfg.topolabel       = topolabel;
    cfg.mixing          = mixing;
    cfg.unmixing        = unmixing;
    cfg.updatesens      = 'no';
    cfg.demean          = 'no';
    cfg.doscale         = 'no';
    cfg.feedback        = 'no';
    
    compdata            = ft_componentanalysis ( cfg, trialdata );
    compdata.grad       = [];
    compdata.elec       = [];
    
    % Corrects the topology maps using the original mixing matrix.
    for mychan = 1: numel ( topolabel )
        
        ftchan = strcmp ( compdata.topolabel, topolabel { mychan } );
        
        mytopo = mixing ( mychan, : );
        compdata.topo ( ftchan, : ) = mytopo;
    end
    
    % Removes the artifact components.
    cfg.component       = cat ( 1, EOGcomp, EKGcomp );
    cleandata           = ft_rejectcomponent ( cfg, compdata, trialdata );
    
    % Gets the list of non-MEG channels (EOG, EKG, etc.).
    MEG                 = find ( ismember ( cleandata.label', ft_channelselection ( { 'MEG' },             cleandata.label ) ) );
    MEGMAG              = find ( ismember ( cleandata.label', ft_channelselection ( { 'MEGMAG' },          cleandata.label ) ) );
    MEGGRAD             = find ( ismember ( cleandata.label', ft_channelselection ( { 'MEGGRAD' },         cleandata.label ) ) );
    EEG                 = find ( ismember ( cleandata.label', ft_channelselection ( { 'EEG' },             cleandata.label ) ) );
    physio              = find ( ismember ( cleandata.label', ft_channelselection ( { 'EOG' 'ECG' 'EMG' }, cleandata.label ) ) );
    
    % Sets the colors according to the channel type.
    chancol             = zeros ( numel ( cleandata.label ), 1 );
    chancol ( MEG )     = 1;
    chancol ( MEGMAG )  = 2;
    chancol ( MEGGRAD ) = 3;
    chancol ( EEG )     = 4;
    chancol ( physio )  = 5;
    
    % Defines the EEG projector.
    projector           = [];
    projector.kind      = 10; % EEG average reference.
    projector.active    = true;
    projector.desc      = 'EEG average reference';
    projector.data      = [];
    projector.data.nrow      = 1;
    projector.data.ncol      = sum ( EEG );
    projector.data.row_names = [];
    projector.data.col_names = trialdata.label ( EEG );
    projector.data.data      = ones ( 1, sum ( EEG ) );
    
    % Displays the MEG data to append or remove artifacts.
    cfg                 = [];
    cfg.channel         = 1: 30;
    cfg.physio          = physio;
    cfg.gradscale       = 0.05;
    cfg.eegscale        = 1e-13;
    cfg.eogscale        = 1e-13;
    cfg.ecgscale        = 1e-13;
    cfg.ylim            = [ -1 1 ] * 4e-12;
    cfg.plotlabels      = 'yes';
    cfg.viewmode        = 'vertical';
    cfg.continous       = 'yes';
    cfg.colorgroups     = chancol;
    cfg.channelcolormap = [   0 114 189;   0 114 189; 162  20  47; 217  83  25; 126  47 142 ] / 255;
    cfg.artfctdef       = artifact;
    cfg.ploteventlabels = 'colorvalue';
    cfg.event           = event;
    
    cfg.badchan         = fileconfig.channel.bad;
    cfg.projinfo        = projector;
    cfg.applyprojector  = true;
    
    cfg                 = my_databrowser ( cfg, cleandata );
    drawnow
     
    % Uses the full list of artifacts.
    artifact        = cfg.artfctdef;
    
    % Gets the previous step artifact definition.
    artinfo         = metadata.artinfo;
    
    % Gets the number of previous and current samples.
    psamples        = 0;
    csamples        = headers.nSamples;
    
    % Corrects the artifact definitions.
    arttypes = fieldnames ( artifact );
    for aindex = 1: numel ( arttypes )
        
        % Gets the list of artifacts of the current type.
        arttype         = arttypes { aindex };
        artifact        = cfg.artfctdef.( arttype ).artifact;
        
        % Resamples and corrects the artifact definitions.
        artifact        = artifact * downrate;
        artifact        = artifact - psamples;
        
        % Removes the artifacts of previous files.
        artifact ( artifact ( :, 2 ) < 1,        : ) = [];
        
        % Removes the artifacts of future files.
        artifact ( artifact ( :, 1 ) > csamples, : ) = [];
        
        % Corrects the first and last artifact, if needed.
        artifact        = max ( artifact, 1 );
        artifact        = min ( artifact, csamples );
        
        % Stores the current artifact definition.
        artinfo.artifact.( arttype ).artifact = artifact;
    end
    
    % Updates the artifacts
    metadata.artinfo ( findex ) = artinfo;
    
    % Update the history
    metadata.history.last_step = 'Artifact revision';
    metadata.history.previous_steps(end+1) = {metadata.history.last_step};
    metadata.history.date = datestr ( now );
    
    
    % Saves the output data.
    save ( '-v6', metadata_file, '-struct', 'metadata' )
end
