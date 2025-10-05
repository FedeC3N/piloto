%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Check the triggers:
% - Are they in conductual and EEG?
% - Are the times the same?
%
% 05/10/2025
% Federico Ramírez-Toraño
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
close all

% Add required folders to PATH
addpath(fullfile('..','shared','fieldtrip-20250928'));
ft_defaults

% Define paths to data
config.path.eeg = fullfile('..','..','data','eeg');
config.path.conductual = fullfile('..','..','data','conductual');

% Find the subjects based on the EEG
files = dir(sprintf('%s/*WM*.cnt',config.path.eeg));
files = {files.name};

for ifile = 1 : numel(files)
    
    % Get the current subject code
    current_subject = files{ifile};
    current_subject = current_subject(1:9);
    fprintf('Working on subject %s\n\n', current_subject)
    
    % Load the triggers in conductual data
    conductual = dir(sprintf('%s/%s_ColorK_MATLAB.mat',config.path.conductual,current_subject));
    if numel(conductual) == 1
        conductual = load(sprintf('%s/%s',config.path.conductual,conductual.name));
    else
        fprintf('  Incorrect number of conductual files: %i\n\n', numel(conductual));
    end
    conductual_triggers = conductual.stim.triggers;
    
    % Load the EEG events and the header
    eeg_dataset = dir(sprintf('%s/%s*WM*.cnt',config.path.eeg,current_subject));
    current_eeg_file = sprintf('%s/%s',eeg_dataset.folder,eeg_dataset.name);
    eeg_triggers = ft_read_event(current_eeg_file);
    header = ft_read_header(current_eeg_file);
    eeg_triggers = convert_triggers_format(eeg_triggers,header);
    
    %%%% Check of conductual trials
    check_conductual_triggers(conductual_triggers);
    
    %%%% Check of eeg trials
    check_eeg_triggers(eeg_triggers,header);
    
    
end


%%%%%%%% FUNCTIONS
function eeg_triggers = convert_triggers_format(eeg_triggers,header)


% Save the original values
dummy = eeg_triggers;

% Create the strcut again
eeg_triggers = [];
eeg_triggers.onset = [dummy.sample] / header.Fs;
eeg_triggers.value = [dummy.value];


end

function check_conductual_triggers(conductual_triggers)

start_experiment_value = 8;
start_experiment_index = conductual_triggers.value == start_experiment_value;
start_onset = conductual_triggers.onset(start_experiment_index);

end_experiment_value = 72;
end_experiment_index = conductual_triggers.value == end_experiment_value;
end_onset = conductual_triggers.onset(end_experiment_index);

start_trial_value = 16;
start_trial_index = conductual_triggers.value == start_trial_value;
num_of_trials = numel(conductual_triggers.value(start_trial_index));

response_one_value = 1;
response_two_value = 2;
response_index = conductual_triggers.value == response_one_value | conductual_triggers.value == response_two_value;
num_of_rseponse = numel(conductual_triggers.value(response_index));

fprintf('%%%%%%%%%%%%%%%%  Check of conductual triggers\n')
fprintf('Num of start triggers (%i): %i\n', start_experiment_value, sum(start_experiment_index))
fprintf('Seconds at start: %.3f\n',start_onset)
fprintf('Minutes at end: %.3f\n',end_onset/60)
fprintf('Num of trials (%i): %i\n', start_trial_value, num_of_trials)
fprintf('Num of responses (1 - 2): %i\n', num_of_rseponse)

fprintf('\n\n\n')


end



function check_eeg_triggers(eeg_triggers,header)

start_experiment_value = 8;
start_experiment_index = eeg_triggers.value == start_experiment_value;
start_onset = eeg_triggers.onset(start_experiment_index);

end_experiment_value = 72;
end_experiment_index = eeg_triggers.value == end_experiment_value;
end_onset = eeg_triggers.onset(end_experiment_index);

start_trial_value = 16;
start_trial_index = eeg_triggers.value == start_trial_value;
num_of_trials = numel(eeg_triggers.value(start_trial_index));

response_one_value = 1;
response_two_value = 2;
response_index = eeg_triggers.value == response_one_value | eeg_triggers.value == response_two_value;
num_of_rseponse = numel(eeg_triggers.value(response_index));

fprintf('%%%%%%%%%%%%%%%%  Check of EEG triggers\n')
fprintf('Num of start triggers (%i): %i\n', start_experiment_value, sum(start_experiment_index))
fprintf('Seconds at start: %.3f\n',start_onset)
fprintf('Minutes at end: %.3f\n',end_onset/60)
fprintf('Num of trials (%i): %i\n', start_trial_value, num_of_trials)
fprintf('Num of responses (1 - 2): %i\n', num_of_rseponse)

fprintf('\n\n\n')


end

