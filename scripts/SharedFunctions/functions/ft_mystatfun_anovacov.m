function [s]= ft_mystatfun_anovacov(cfg,data,design)

% Data is a nfuentes x n subjects vector with the data to analyze

% cfg.grupo is a 1x2 cell where:
% (1,1 is an empty vector wich will be filled with design (independent
% variable)
% (1,2) is 1xnsubject vector with the cov1
% (1,3) is 1xnsubject vector with the cov2 QUITA

% design is 1xnsubject vector with the independent variable
if ~isfield(cfg,'group')
    df1 = length(unique(design))-1;
    df2 = length(design)-1;
    %Critval calculado para cfg.tail = 0 (suponiendo doble cola)
    s.critval = [finv(1-cfg.alpha, df1, df2)];
else
    df1 = length(unique(design))-1;
    df2 = length(design)-1;
    %Critval calculado para cfg.tail = 0 (suponiendo doble cola)
    s.critval = [finv(1-cfg.alpha, df1, df2)];
    
    grupo = cfg.group;
    grupo{1}=design';
    
    
    [ ~, table ] = my_anovan ( data', grupo, 'continuous', 2, 'display', 'off' );
%     [~, table] = my_anovan ( data', grupo (1), 'display', 'off' );
    
    s.stat = squeeze(table(2,6,:));
    s.prob = squeeze(table(2,7,:));
%         s.stat = squeeze(table(3,6,:));
%     s.prob = squeeze(table(3,7,:));
    
    s.stat = [s.stat{:}]';
    s.prob = [s.prob{:}]';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EDAD
%     s.stat_edad = squeeze(table(3,6,:));
%     s.prob_edad = squeeze(table(3,7,:));
%     
%     s.stat_edad = [s.stat_edad{:}]';
%     s.prob_edad = [s.prob_edad{:}]';
    
    % s.prob lo mismo
    % for ifu = 1:size(data,1)%28968
    % aux = data(ifu,:);
    % [p,table] = anovan(aux, grupo, 'continuous', [2], 'display','off');
    % s.stat(ifu,1) = table{2,6};
    % s.prob(ifu,1) = table{2,7};
    % end
end
end

