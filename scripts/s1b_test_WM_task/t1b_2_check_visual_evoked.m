%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Check visual evoked
%
% 05/10/2025
% Federico Ramirez-Torano
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
close all

% Add required folders to PATH
addpath(fullfile('..','shared','fieldtrip-20250928'));
ft_defaults

% Define paths to data
config.path.eeg = fullfile('..','..','data','eeg','raw');
config.path.conductual = fullfile('..','..','data','conductual');

% Find the subjects based on the EEG
files = dir(sprintf('%s/*WM*.cnt',config.path.eeg));
% files = dir(sprintf('%s/*prueba_64ch*.cnt',config.path.eeg));
files = {files.name};

for ifile = 3 : numel(files)
    
    % Get the current subject code
    current_subject = files{ifile};
    current_subject = current_subject(1:9);
%     current_subject = current_subject(1:19);
    fprintf('Working on subject %s\n\n', current_subject)
    
    % Load the triggers in conductual data
    conductual = dir(sprintf('%s/%s_ColorK_MATLAB.mat',config.path.conductual,current_subject));
    if numel(conductual) == 1
        conductual = load(sprintf('%s/%s',config.path.conductual,conductual.name));
    else
        fprintf('  Incorrect number of conductual files: %i\n\n', numel(conductual));
    end
    conductual_triggers = conductual.stim.triggers;
    
    % Load the EEG, events and the header
    eeg_dataset = dir(sprintf('%s/%s',config.path.eeg,files{ifile}));
    current_eeg_file = sprintf('%s/%s',eeg_dataset.folder,eeg_dataset.name);
    eeg_triggers = ft_read_event(current_eeg_file);
    header = ft_read_header(current_eeg_file);
    
    % Align the triggers
    t_diff = eeg_triggers(end).sample/header.Fs - conductual_triggers.onset(end);
    conductual_triggers.onset = conductual_triggers.onset + t_diff;
    
    % trigger index
    trigger_of_interest = [16, 24, 32, 40, 48];
    trigger_of_interest = [40,48];
    trigger_index = ismember(conductual_triggers.value,trigger_of_interest);
    onset = round(conductual_triggers.onset(trigger_index)*header.Fs);
    values = conductual_triggers.value(trigger_index);
    
    % Define the trial
    prestim = 0.200; %seconds
    poststim = 0.800; %seconds
    trlbeg = onset - round(prestim*header.Fs);
    trlend = onset + round(poststim*header.Fs);
    offset = -round(prestim*header.Fs);
    trl = [trlbeg(:), trlend(:) repmat(offset,numel(trlbeg),1), values'];
    
    % Read the EEG with this trl definition
    cfg = [];
    cfg.dataset = current_eeg_file;
    cfg.channel = 'EEG';
    cfg.trl = trl;
    cfg.hpfilter = 'yes'; cfg.hpfreq = 2;
    cfg.lpfilter = 'yes'; cfg.lpfreq = 45;
    cfg.reref = 'yes'; cfg.refchannel = 'all';
    data = ft_preprocessing(cfg);
    
    % clean data
    cfg = [];
    cfg.layout = 'natmeg_customized_eeg1005.lay';
    cleaned_data_ERP = ft_rejectvisual(cfg,data);
    
    % Estimate ERP
    cfg = [];
    ERP = ft_timelockanalysis(cfg,cleaned_data_ERP);
    
    % Create the layout
    templ = fullfile('..','shared','fieldtrip-20250928','template','electrode','standard_1020.elc');
    elec = ft_read_sens(templ);
    [common, ia, ib] = intersect(elec.label,cleaned_data_ERP.label,'stable');
    elec.label = elec.label(ia);
    elec.chanpos = elec.chanpos(ia,:);
    elec.elecpos = elec.elecpos(ia,:);
    layout = ft_prepare_layout(struct('elec',elec));
    
    
    % Plot
    figure
    cfg = [];
    cfg.layout = layout;
    ft_multiplotER(cfg,ERP);
    
    
    
%     % Plot 
%     cfg = [];
%     cfg.trl = trl;
%     cfg.channel = data.label(1:10);
%     cfg.continuous = 'no';
%     cfg.viewmode = 'vertical';
%     ft_databrowser(cfg,data);



    
    
    
end

    
function data_clena = quick_automatic_clean(data)

% 1) Flatline / insanely noisy channels (auto flag, then repair)
%    We compute per-channel robust SD; flag extreme z-scores.
tmp = ft_timelockanalysis(struct('keeptrials','yes'), data);
dat = tmp.trial;  % trials x chans x time
chanSD = squeeze(std(dat,0,[1 3]));             % SD across trials & time
m = median(chanSD);  s = 1.4826*median(abs(chanSD-m));   % robust MAD->SD
z = (chanSD - m)/s;
badChan = find(abs(z) > 4 | chanSD < 1e-7);     % threshold can be tightened/loosened

% 2) Muscle artifacts (110–140 Hz z-threshold on envelope)
cfg = [];
cfg.trl       = data.sampleinfo;  % in case ft_artifact_* needs it
cfg.continuous= 'no';
cfg.artfctdef.muscle.channel  = 'EEG';
cfg.artfctdef.muscle.bpfilter = 'yes';
cfg.artfctdef.muscle.bpfreq   = [80 140];
cfg.artfctdef.muscle.hilbert  = 'yes';
cfg.artfctdef.muscle.cutoff   = 10;         % z-value cutoff
[cfg, art_muscle] = ft_artifact_muscle(cfg, data);

% 3) Blink/EOG with zvalue (use EOG; if missing, fall back to frontal proxies)
allLabels = data.label;
VEOGuse = intersect({'Fp2','Fpz','Fp1'}, allLabels, 'stable');
HEOGuse = intersect({'F8','F7','AF7','AF8'}, allLabels, 'stable');

cfg = [];
cfg.trl       = data.sampleinfo;
cfg.continuous= 'no';
cfg.artfctdef.zvalue.channel = unique([VEOGuse(:); HEOGuse(:)])';
cfg.artfctdef.zvalue.cutoff  = 10;          % z cutoff for blinks/saccades
cfg.artfctdef.zvalue.demean  = 'yes';
cfg.artfctdef.zvalue.hilbert = 'no';
cfg.artfctdef.zvalue.lpfilter= 'yes';
cfg.artfctdef.zvalue.lpfreq  = 15;         % blink energy <~15 Hz
[cfg, art_eog] = ft_artifact_zvalue(cfg, data);

% 4) Jumps (step-like glitches)
cfg = [];
cfg.trl       = data.sampleinfo;
cfg.continuous= 'no';
cfg.artfctdef.jump.channel = 'EEG';
cfg.artfctdef.jump.medianfilter = 'yes';
cfg.artfctdef.jump.medianfiltord= 9;
cfg.artfctdef.jump.cutoff       = 10;       % z cutoff
[cfg, art_jump] = ft_artifact_jump(cfg, data);

% Combine artifact segments and reject entire trials that overlap artifacts
cfg = [];
cfg.artfctdef.reject = 'complete';   % drop trials that contain artifacts
cfg.artfctdef.muscle.artifact = art_muscle;
cfg.artfctdef.eog.artifact    = art_eog;
cfg.artfctdef.jump.artifact   = art_jump;
data_clean = ft_rejectartifact(cfg, data);

% Report rejections
nAll  = numel(data.trial);
nKeep = numel(data_clean.trial);
fprintf('Trials kept: %d / %d (%.1f%%)\n', nKeep, nAll, 100*nKeep/nAll);


end
