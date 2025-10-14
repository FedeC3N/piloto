clc
clear
close all

% Adds the functions folders to the path.
addpath ( '../SharedFunctions/fieldtrip-20220626/');
addpath ( '../SharedFunctions/functions/');
addpath ( '../SharedFunctions/mne_silent/');
addpath ( '../SharedFunctions/functions_eep/');
ft_defaults

% Select user
fid = fopen('./user.txt','r');
user = fgetl(fid);
fclose(fid);

% Paths
config.path.raw  = '../../databases/AI_Mind_database';
config.path.derivatives = sprintf('../../databases/AI_Mind_database/derivatives/%s',user);

% Action when the task have already been processed.
config.overwrite      = false;

% Sets the segmentation options.
config.trialfun       = 'restingSegmentation';
config.segment        = 4.0;
config.padding        = 2.0;
config.overlap        = 0.0;

% Sets the sketch building options.
config.freqband       = [ 2 45 ];

config.channel.groups = { 'EEG' };
config.channel.ignore = {};

% Lists the files.
dataset_file = sprintf('%s/dataset_%s.mat',config.path.derivatives,user);
dataset    = struct2array ( load ( dataset_file ) );

% Goes through each subject and task.
for sindex = 1: numel ( dataset )
    
    % Read the metadata
    current_folder = dataset(sindex).folder;
    current_file = dataset(sindex).file;
    current_subject = dataset(sindex).subject;
    current_session = dataset(sindex).session;
    current_task = dataset(sindex).task;
    fprintf('Current file sub-%s_ses-%s_task-%s\n', current_subject,current_session,current_task);
    
    % Check if exist metadata file (mandatory)
    current_path_trl = sprintf('%s/sub-%s/ses-%s/eeg',...
            config.path.derivatives,current_subject,current_session);
    metadata_file   =  sprintf ( '%s/sub-%s_ses-%s_task-%s_desc-trl_eeg.mat', current_path_trl, ...
        current_subject, current_session, current_task );
    if ~exist(metadata_file)
        fprintf('Metadata file not created. Skip \n\n');
    end
    
    % Load metadata
    metadata = load (metadata_file);
    dummy_path_current_file = sprintf('%s/%s',current_folder,current_file);
    
    % Checks that the selected channel group is present in the SOBI data.
    if ~any ( ismember ( config.channel.groups, fieldnames ( metadata.compinfo.SOBI ) ) )
        fprintf ( 1, 'Ignoring %s (no SOBI information for the selected channel groups).\n', msgtext );
        continue
    end
    
    % Check if overwrite
    outname             = sprintf ( '%s/sub-%s_ses-%s_task-%s_desc-sketch_eeg.mat', current_path_trl, metadata.subject, metadata.session, metadata.task ); 
    if exist ( outname, 'file' ) && ~config.overwrite
        fprintf ( 1, '  Already calculated).\n');
        continue
    end
    
    % Gets the files length as a vector.
    headers             = [ metadata.fileinfo.header ];
    samples             = [ headers.nSamples ];
    offsets             = cat ( 2, 0, cumsum ( samples) );
    
    % Reserves memory for the trialdata.
    trialdefs           = cell ( numel ( metadata.fileinfo ), 1 );
    trialfiles          = cell ( numel ( metadata.fileinfo ), 1 );
    trialpads           = cell ( numel ( metadata.fileinfo ), 1 );
    trialtypes          = cell ( numel ( metadata.fileinfo ), 1 );
    erfdatas            = cell ( numel ( metadata.fileinfo ), 1 );
    freqdatas           = cell ( numel ( metadata.fileinfo ), 1 );
    
    % Goes through all the data files.
    for findex = 1: numel ( metadata.fileinfo )
        
        % Loads the current epochs and artifact definitions.
        fileinfo             = metadata.fileinfo ( findex );
        artinfo              = metadata.artinfo  ( findex );
        dummy_path_current_file = sprintf('%s/%s',current_folder,fileinfo.file);

        
        if ~exist ( dummy_path_current_file, 'file' )
            fprintf ( 1, '  Ignoring trial definition for file %i (%s file not found).\n', findex, basename );
            continue
        end
        
        
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
        
        
        % Marks the artifacted trials.
        
        
        % Writes the trial's metadata.
        trialtype             = zeros ( numel ( erfdata.trial ), 1 );
        trialpad              = zeros ( numel ( erfdata.trial ), 1 );
        trialfile             = repmat ( findex,   numel ( erfdata.trial ), 1 );
        
        
        % Stores the epoch data.
        trialdefs  { findex } = trialdef;
        trialfiles { findex } = trialfile;
        trialpads  { findex } = trialpad;
        trialtypes { findex } = trialtype;
        erfdatas   { findex } = erfdata;
        freqdatas  { findex } = freqdata;
        
        clear erfdata
        clear freqdata
    end
    
    % Removes the empty files.
    empty                 = cellfun ( @isempty, erfdatas );
    
    erfdatas   ( empty )  = [];
    freqdatas  ( empty )  = [];
    
    trialdefs  ( empty )  = [];
    trialfiles ( empty )  = [];
    trialpads  ( empty )  = [];
    trialtypes ( empty )  = [];
    
    trialdef              = cat ( 1, trialdefs  {:} );
    trialfile             = cat ( 1, trialfiles {:} );
    trialpad              = cat ( 1, trialpads  {:} );
    trialtype             = cat ( 1, trialtypes {:} );
    
    clear trialdefs
    clear trialfiles
    clear trialpads
    clear trialtypes
    
    
    if isempty ( erfdatas )
        fprintf ( 1, '  No data for subject ''%s'' task ''%s''. Skipping.\n', metadata.subject, metadata.task );
        continue
    end
    
    fprintf ( 1, '  Merging all the data files.\n' );
    
    % Merges the epoch data.
    if numel ( erfdatas ) > 1
        
        % Merges the ERF.
        cfg           = [];
        
        erfdata       = myft_appenddata ( cfg, erfdatas  {:} );
        
        % Merges the frequency data.
        cfg.parameter = 'fourierspctrm';
        
        freqdata      = ft_appendfreq ( cfg, freqdatas {:} );
        
        % Removes the 'cfg' field.
        freqdata      = rmfield ( freqdata,  'cfg' );
    else
        
        erfdata       = erfdatas   {1};
        freqdata      = freqdatas  {1};
    end
    clear erfdatas
    clear freqdatas
    
    
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
        outname             = sprintf ( '%s/sub-%s_ses-%s_task-%s_desc-sketch_eeg.mat', current_path_trl, metadata.subject, metadata.session, metadata.task );
        
        if exist ( outname, 'file' ) && ~config.overwrite
            fprintf ( 1, '  Ignoring channel group ''%s'' (already calculated).\n', channel );
            continue
        end
        
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
        groupinfo.subject     = metadata.subject;
        groupinfo.session     = metadata.session;
        groupinfo.task        = metadata.task;
        groupinfo.stage       = metadata.stage;
        groupinfo.channel     = channel;
        groupinfo.fileinfo    = metadata.fileinfo;
        groupinfo.chaninfo    = metadata.chaninfo;
        groupinfo.artinfo     = metadata.artinfo;
        groupinfo.compinfo    = compinfo;
        groupinfo.trialinfo   = trialinfo;
        groupinfo.erfdata     = grouperf;
        groupinfo.freqdata    = groupfreq;
        groupinfo.cleaninfo   = cleaninfo;
        
        % Saves the current group epoch data.
        save ( '-v6', outname, '-struct', 'groupinfo' );
    end
end
