%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Quick numbers for beh pilot
%
% 01/10/2025
% Federico Ramirez-Torano
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
close all

% Paths
config.path.beh  = fullfile('..','..','data','beh');

% Config
index_load = [4, 6, 8];
num_blocks = 5;
colors = [0    0.4470    0.7410;...
    0.8500    0.3250    0.0980;...
    0.9290    0.6940    0.1250;...
    0.4940    0.1840    0.5560;...
    0.4660    0.6740    0.1880;...
    0.3010    0.7450    0.9330;...
    0.6350    0.0780    0.1840];

% Get the subjects' file
files = dir(fullfile(config.path.beh,'HIQ_0*_ColorK_MATLAB.mat'));

% Estimate the beh information of interest
block_accuracy = nan(numel(index_load),num_blocks,numel(files));
K = nan(numel(index_load),numel(files));
for ifile = 1 : numel (files)
    
    % Load the current file
    current_file = fullfile(files(ifile).folder,files(ifile).name);
    load(current_file);
    
    fprintf('WORKING ON %s \n',files(ifile).name);
    
    % Get the index of the load
    accuracy = stim.accuracy;
    task_load = stim.setSize;
    
    % Global accuracy
    global_accuracy = [...
        sum(accuracy(task_load == index_load(1)))/ numel(accuracy(task_load == index_load(1))),...
        sum(accuracy(task_load == index_load(2)))/ numel(accuracy(task_load == index_load(2))),...
        sum(accuracy(task_load == index_load(3)))/ numel(accuracy(task_load == index_load(3)))];
    
    fprintf(1,'Global Index 4 accuracy: %.3f\n', global_accuracy(1))
    fprintf(1,'Global Index 6 accuracy: %.3f\n', global_accuracy(2))
    fprintf(1,'Global Index 8 accuracy: %.3f\n\n', global_accuracy(3))
    
    % Accuracy per block
    for iload = 1 : numel(index_load)
        current_block_accuracy = accuracy(task_load == index_load(iload));
        current_block_accuracy = reshape(current_block_accuracy,[],5);
        block_accuracy(iload,:,ifile) = sum(current_block_accuracy,1)/size(current_block_accuracy,1);
        %         plot(ax1,1:5,block_accuracy(iload,:),'*-')
    end
    %     legend('Set 4', 'Set 6', 'Set 8')
    %     ylim([0 1])
    %     xticks([1 2 3 4])
    %     xlabel('Block number')
    %     ylabel('Accuracy')
    
    % K coefficient
    [current_K,setSize] = computeK(stim);
    K(:,ifile) = current_K;
    
    % Plot
    %     plot(setSize, K, '-o', 'LineWidth', 2, 'MarkerSize', 8);
    %     ylim([0 max(K)*1.2]);
    %     xlim([min(setSize)-0.5, max(setSize)+0.5]);
    %
    %     xlabel('Set Size (# of items)');
    %     ylabel('Estimated K');
    %     title('Visual Working Memory Capacity (K) by Set Size');
    %     grid on;
    %
    %     % Estética tipo paper
    %     set(gca, 'FontSize', 12);
    
    
    % Time spent on the task
    fprintf(1,'Time elapsed: %.3f minutes. \n\n\n\n', stim.triggers.onset(end)/60)
    
end


% Plot the information of interest
% Plot global accuracy
global_accuracy = squeeze(mean(block_accuracy,2));
figure
for ifile = 1 : size(global_accuracy,2)
    
    plot(index_load,global_accuracy(:,ifile),'*-')
    hold on
    ylim([0 1])
    xticks(index_load)
    xlabel('Load')
    ylabel('Accuracy')
    
end

% Plot K
figure
for ifile = 1 : size(global_accuracy,2)
    
    plot(index_load,K(:,ifile),'*-')
    hold on
    xticks(index_load)
    xlabel('Load')
    ylabel('K')
    
end



% Plot accuracy per block
% for iload = 1 : numel(index_load)
%     figure
%     current_y = squeeze(block_accuracy(iload,:,:));
%     plot(1:num_blocks,current_y,'*-','Color',colors(iload,:))
%     hold on
%     ylim([0 1])
%     xticks([1 2 3 4 5])
%     xlabel('Block number')
%     ylabel('Accuracy')
%     title(sprintf('Load %i', index_load(iload)))
%     
% end





%%% LOCAL FUNCTION
function [K, uniqueSetSizes] = computeK(stim)
% computeK_bySetSize - Calcula la capacidad de memoria visual (K) agrupada por tamaño del set
%
% Input:
%   stim.setSizes : matriz [120 x N] con los tamaños del set por ensayo
%   stim.accuracy : matriz [120 x N] con 1 = respuesta correcta, 0 = incorrecta
%   stim.change   : matriz [120 x N] con 1 = hubo cambio, 0 = no hubo cambio
%
% Output:
%   K               : vector con valores de K para cada tamaño de set único
%   uniqueSetSizes  : vector con los tamaños de set correspondientes a cada K

% Aplanar todas las matrices en vectores
setSize = stim.setSize(:);
accuracy = stim.accuracy(:);
change = stim.change(:);

% Identificar los tamaños únicos de set
uniqueSetSizes = unique(setSize);
K = zeros(size(uniqueSetSizes));

for i = 1:length(uniqueSetSizes)
    S = uniqueSetSizes(i);
    
    % Filtrar ensayos que tienen este set size
    idx = setSize == S;
    
    acc = accuracy(idx);
    chg = change(idx);
    
    % Tasa de aciertos
    H = mean(acc(chg == 1));
    
    % Tasa de falsas alarmas
    F = 1 - mean(acc(chg == 0));
    
    % Calcular K
    K(i) = S * (H - F);
end
end


