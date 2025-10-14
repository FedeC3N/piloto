function [s]= ft_mystatfun_anovacov_2_factors(cfg,data,design)

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
    
    grupo = cfg.group;
    %     grupo{1}=design;
    
    if isempty(cfg.permute)
        if isfield(cfg,'model')
            [~, table] = my_anovan ( data', grupo,...
                'model',cfg.model,'display', 'off' );
        else
            [~, table] = my_anovan ( data', grupo ,...
                'display', 'off' );
        end
        
        s.stat = nan(size(table,3),size(cfg.model,1));
        s.prob = nan(size(table,3),size(cfg.model,1));
        s.critval = nan(1,size(cfg.model,1));
            
        
        % FACTOR 1
        dummy_stat = squeeze(table(2,6,:));
        dummy_prob = squeeze(table(2,7,:));
        s.stat(:,1) = [dummy_stat{:}]';
        s.prob(:,1) = [dummy_prob{:}]';
        s.critval(1) = finv(1-cfg.alpha,table{2,3,1},table{5,3,1});
        clear dummy_stat dummy_prob 
        
        % FACTOR 2
        dummy_stat = squeeze(table(3,6,:));
        dummy_prob = squeeze(table(3,7,:));
        s.stat(:,2) = [dummy_stat{:}]';
        s.prob(:,2) = [dummy_prob{:}]';
        s.critval(2) = finv(1-cfg.alpha,table{3,3,1},table{5,3,1});
        clear dummy_stat dummy_prob
        
        % INTERACCION
        dummy_stat = squeeze(table(4,6,:));
        dummy_prob = squeeze(table(4,7,:));
        s.stat(:,3) = [dummy_stat{:}]';
        s.prob(:,3) = [dummy_prob{:}]';
        s.critval(3) = finv(1-cfg.alpha,table{4,3,1},table{5,3,1});
        clear dummy_stat dummy_prob
        
        
    else
        
        s.stat = nan(size(data,1),size(cfg.model,1));
        s.prob = nan(size(data,1),size(cfg.model,1));
        s.critval = nan(1,size(cfg.model,1));
        
        % Permutate FAM vs CN
        newg1 = grupo{1};
        newg1 = newg1(randperm(numel(newg1)));        
        [~,table,~] = my_anovan(data',{newg1 grupo{2}}...
            , 'model',cfg.model,...
            'display','off');
        dummy_stat = squeeze(table(2,6,:));
        dummy_prob = squeeze(table(2,7,:));
        s.stat(:,1) = [dummy_stat{:}]';
        s.prob(:,1) = [dummy_prob{:}]';
        s.critval(1) = finv(1-cfg.alpha,table{2,3,1},table{5,3,1});
        clear newg1 table dummy_stat dummy_prob
        
        % Permutate APOE+ vs APOE-
        newg2 = grupo{2};
        newg2 = newg2(randperm(numel(newg2)));        
        [~,table,~] = my_anovan(data',{grupo{1} newg2},...
            'model',cfg.model,...
            'display','off');
        dummy_stat = squeeze(table(3,6,:));
        dummy_prob = squeeze(table(3,7,:));
        s.stat(:,2) = [dummy_stat{:}]';
        s.prob(:,2) = [dummy_prob{:}]';
        s.critval(2) = finv(1-cfg.alpha,table{3,3,1},table{5,3,1});
        clear newg2 table dummy_stat dummy_prob
        
        % Permutate the interaction
        newg1 = grupo{1};
        newg2 = grupo{2};
        dummy_index_perm = randperm(numel(newg1));
        newg1 = newg1(dummy_index_perm);
        newg2 = newg2(dummy_index_perm);
        [~,table,~] = my_anovan(data',{newg1 newg2},...
            'model',cfg.model,...
            'display','off');
        dummy_stat = squeeze(table(4,6,:));
        dummy_prob = squeeze(table(4,7,:));
        s.stat(:,3) = [dummy_stat{:}]';
        s.prob(:,3) = [dummy_prob{:}]';
        s.critval(3) = finv(1-cfg.alpha,table{4,3,1},table{5,3,1});
        clear newg1 newg2 table dummy_stat dummy_prob dummy_index_perm
        
        
        
    end
end
end

