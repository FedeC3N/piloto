%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Quick numbers for pilot
%
% 01/10/2025
% Federico Ramírez-Toraño
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
close all

load('../../data/conductual/HIQ_022_0_ColorK_MATLAB.mat');

% Get the index of the load
accuracy = stim.accuracy;
task_load = stim.setSize;
index_load = [4, 6, 8];

% Global accuracy
global_accuracy = [...
sum(accuracy(task_load == index_load(1)))/ numel(accuracy(task_load == index_load(1))),...
sum(accuracy(task_load == index_load(2)))/ numel(accuracy(task_load == index_load(2))),...
sum(accuracy(task_load == index_load(3)))/ numel(accuracy(task_load == index_load(3)))];

fprintf(1,'Global Index 4 accuracy: %.3f\n', global_accuracy(1))
fprintf(1,'Global Index 6 accuracy: %.3f\n', global_accuracy(2))
fprintf(1,'Global Index 8 accuracy: %.3f\n\n', global_accuracy(3))

% Accuracy per block
figure
block_accuracy = nan(3,5);
for iload = 1 : 3
  current_block_accuracy = accuracy(task_load == index_load(iload));
  current_block_accuracy = reshape(current_block_accuracy,[],5);
  block_accuracy(iload,:) = sum(current_block_accuracy,1)/size(current_block_accuracy,1);
  plot(1:5,block_accuracy(iload,:),'*-')
  hold on
end
legend('Set 4', 'Set 6', 'Set 8')
ylim([0 1])
xticks([1 2 3 4])
xlabel('Block number')
ylabel('Accuracy')

% K coefficient
[K,setSize] = computeK(stim);
figure;
plot(setSize, K, '-o', 'LineWidth', 2, 'MarkerSize', 8);
ylim([0 max(K)*1.2]);
xlim([min(setSize)-0.5, max(setSize)+0.5]);

xlabel('Set Size (# of items)');
ylabel('Estimated K');
title('Visual Working Memory Capacity (K) by Set Size');
grid on;

% Estética tipo paper
set(gca, 'FontSize', 12);


% Time spent on the task
fprintf(1,'Time elapsed: %.3f minutes. \n\n', stim.triggers.onset(end)/60)



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


