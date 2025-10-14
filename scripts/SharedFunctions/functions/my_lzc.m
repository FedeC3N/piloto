function LZC = my_lzc(data)

% Lempel-Ziv complexity of a given time data


% Data dimensions
num_channels = size(data,1);
num_timepoints = size(data,2);
num_trials = size(data,3);



% We select a threshold to do the coarse-graining of the original time
% data in a sequence of 2 symbols. We choose the median due to its
% robustness to outliers (one could choose the mean too)
% Since I'm using the same dictionary for all the trials of a source, I use
% the median of along ALL the trials
dummy = data;
dummy = reshape(dummy,size(dummy,1),[]);
threshold = median(dummy,2);
threshold = repmat(threshold,1,num_timepoints,num_trials);

% LZC upper bound
b = size(dummy,2)/log2(size(dummy,2));

% Transform the data into symbols
data = int8(data >= threshold);

% We set some initial values to compute LZ complexity
c = ones(num_channels,1);	% Initial complexity

fprintf('Channel ')
for ichannel = 1:num_channels
    
    
    if ichannel>1
        for j=0:log10(ichannel-1)
            fprintf('\b'); % delete previous counter display
        end
    end
    fprintf('%d', ichannel);    
    
    S = data(ichannel,1,1);
    
    for itrial = 1:num_trials
        
        current_data = data(ichannel,:,itrial);
        
%         S = current_data(1);
        Q = current_data(2);
        
        for itime = 2:num_timepoints
            
            SQ = cat(2,S,Q);
            SQ_pi = SQ(1:end-1);
            
            SQ_pi_str = string(SQ_pi);
            SQ_pi_str = join(SQ_pi_str,2);
            Q_str = string(Q);
            Q_str = join(Q_str,2);
            
            if ~contains(SQ_pi_str,Q_str)
                c(ichannel) = c(ichannel)+1;				% We update the complexity values, as we have found a new subsequence
                if (itime+1)>num_timepoints
                    break;
                else
                    S = cat(2,S,Q);		% We build S
                    Q = current_data(itime+1);		% We update Q
                end
            else
                if (itime+1)>num_timepoints
                    break;
                else
                    Q = cat(2,Q,current_data(itime+1));	% We update Q
                end
            end
        end        
        
        
    end
    
    
    
end

LZC = c/b;

end