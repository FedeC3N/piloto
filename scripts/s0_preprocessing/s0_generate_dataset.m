%{

Generate a dataset to work with.
Generate the metadata for the first time for each recording.

@author: Fede
03/09/2025

%}

clear
clc
restoredefaultpath

%%% 
% CONFIG
% Adds the functions folders to PATH.
addpath(fullfile('..','SharedFunctions','functions'));
addpath(fullfile('..','SharedFunctions','mne_silent'));
addpath(fullfile('..','SharedFunctions','functions_eep'));
addpath(fullfile('..','SharedFunctions','fieldtrip-20250928'));
addpath(fullfile('..','SharedFunctions','fieldtrip-20250928','external','eeprobe'))
ft_defaults
clc

% Paths
config.path.raw  = fullfile('..','..','data','eeg','raw');
config.path.metadata = fullfile('..','..','data','eeg','metadata');
config.path.patt = '*.cnt';

% Action when the task has already been processed.
config.overwrite = false;

% Defines the physiological channels (later will be relabeled).
config.physio.EOG     = { 'VEOGL' };
config.physio.EKG     = {};
config.physio.EMG     = {};

%%%%

% Create folders if they don't exist
if ~exist(config.path.metadata,'dir')
    mkdir(config.path.metadata)
end

% Lists the files in the EEG folder.
files = dir (sprintf('%s/%s',config.path.raw,config.path.patt) );
files = {files.name};

% If there is a previous dataset, load the dataset to update it
outfile_dataset  = fullfile('..','..','data','eeg','dataset.mat');
if exist(outfile_dataset,'file')
    load(outfile_dataset);
    old_dataset = dataset; 
    clear dataset;
    already_files = {old_dataset.file};
else
    already_files = {};
end
dataset    = struct ('folder', [], 'file', [], 'subject', [], 'session',[],'task', []);
for ifile = 1 : numel(files)
    
    % print()
    current_folder = config.path.raw;
    current_file = files {ifile};
    fprintf(1,'Working on file %s\n', current_file);
    
    % If the current file is already in the dataset, copy the original
    if ismember(current_file,already_files) && ~config.overwrite

        index = find(ismember(already_files,current_file));
        dataset(ifile) = old_dataset(index);
        continue
        
    end
    
    % Extracts the information.
    pattern = '^(HIQ_\d+)_([0-9]+)_([A-Z0-9]+)_';
    tokens = regexp(current_file, pattern, 'tokens');
    if ~isempty(tokens)
        current_subject = tokens{1}{1}; 
        current_session = tokens{1}{2};
        current_task = tokens{1}{3}; 
    end
    
    % Save the information in the dataset
    dataset(ifile).folder = current_folder;
    dataset(ifile).file = current_file;
    dataset(ifile).subject = current_subject;
    dataset(ifile).session = current_session;
    dataset(ifile).task = current_task;
    
    % Output metadata information
    fileinfo = [];
    artinfo = [];
    chaninfo = [];
    history = [];
    
    % fileinfo
    % 'file','subject','task','stage','begtime','endtime','header','event'
    fileinfo.file = current_file;
    fileinfo.subject = current_subject;
    fileinfo.session = current_session;
    fileinfo.task = current_task;
    
    % Read the header
    dummy_path_current_file = sprintf('%s/%s',current_folder,current_file);
    header = ft_read_header(dummy_path_current_file);
    fileinfo.header = header;
    
    % Look for events
    event = ft_read_event(dummy_path_current_file, 'header',header);
    fileinfo.event = event;
    
    % There is no trigger for begtime or endtime
    fileinfo.begtime = floor(1/header.Fs);
    fileinfo.endtime = floor(header.nSamples/header.Fs);
    
    % artifact
    % 'eog','muscle','jump','visual'
    artinfo.artifact.eog.artifact = zeros(0,2);
    artinfo.artifact.muscle.artifact = zeros(0,2);
    artinfo.artifact.jump.artifact = zeros(0,2);
    artinfo.artifact.visual.artifact = zeros(0,2);
    
    % chaninfo
    chaninfo.bad      = {''};
    
    % history
    history.last_step         = 'Creating Metadata';
    history.date         = datestr ( now );
    history.previous_steps      = {history.last_step};
    
    % Sets the task information.
    metadata          = [];
    metadata.subject  = current_subject;
    metadata.session  = current_session;
    metadata.task     = current_task;
    metadata.fileinfo = fileinfo;
    metadata.chaninfo = chaninfo;
    metadata.history  = history;
    metadata.artinfo = artinfo;
    
    % Saves the output data.
    filename = sprintf ('%s_%s_%s_metadata.mat',  ...
        current_subject, current_session, current_task );
    outfile   = fullfile(config.path.metadata,filename);
    save ( '-v6', outfile, '-struct', 'metadata' );
        
    
    
    
end

% Saves the dataset structure.
save ( '-v6', outfile_dataset, 'dataset' )

