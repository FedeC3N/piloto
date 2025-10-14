function [peaks, peakfittedok]=search_1peak(foi,data,rangepeak,options)

% function that fits a peak with power law background to the powerspectrum
% fits according to log(pow)= B - C*log(f) + A*exp(-((f-fp)/sp).^2)
% and finds the corresponding parameters B, C, A, fp and sp

% INPUTS
%   foi= 1x nfoi frequencies of interest
%   data= ntrials x nfoi = powerspectrum for each foi and trial. The peak
%      can also be fitted to a single or average trial, but in that case 
%      the intertrial variability cannot be estimated and the output 
%      peakfittedok is meaningless
%   rangepeak = [fmin fmax] a-priori range (in Hz) for the peak (ex: [5-14])


% OUTPUTS
%   peaks: structure containing the values of the fitted parameters B, C, 
%       A, fp and sp
%   peakfittedok: true or false, depending on whether the amplitude of the
%       peak exceeds the rms intertrial variability


% implemented according to recommendations in:
%    Lodder, S. S., & van Putten, M. J. (2011). Automated EEG analysis:
%        Characterizing the posterior dominant rhythm. Journal of
%        neuroscience methods, 200(1), 86-93.
%    Chiang, A. K. I., Rennie, C. J., Robinson, P. A., Roberts, J. A., 
%        Rigozzi, M. K., Whitehouse, R. W., ... & Gordon, E. (2008). 
%        Automated characterization of multiple alpha peaks in multi-site
%        electroencephalograms. Journal of neuroscience methods, 168(2),
%        396-411.


peaks=[];
logpf=log(trimmean(data,10,'round',1)); %Average over trials --> logpf: 1 x ntrials
Ntrials=size(data,1);

idrangepeak=find(and(foi>rangepeak(1),foi<rangepeak(2))); %range for the peak frequency


% 1 - Estimates background spectrum logpbackground = B - C*log(f)
pfminusbackground = @(p)(logpf + p(1)*log(foi) - p(2)  );
minB=polyfit(foi,logpf,0);
maxB=polyfit(foi,logpf+2*log(foi),0);
[p,~] = lsqnonlin(pfminusbackground,[0 mean(logpf)],[0 minB],[2 maxB],options);
peaks.C=p(1); peaks.B=p(2);

% 2 - Fits peak

% 2a-initial value for peak frequency and amplitude in the optimization algorithm
[vmax,imax]=max(logpf(idrangepeak)-(peaks.B - peaks.C*log(foi(idrangepeak))));
f1=foi(idrangepeak(imax));

% 2b - finds peak parameters with lsqnonlin optimization
% C stays fixed from the previous optimization
myerror = @(x)(abs(x(1) - peaks.C*log(foi) + x(2)*exp(-((foi-x(4))/x(3)).^2) - logpf));
[x,~] = lsqnonlin(myerror,[peaks.B; vmax; 1; f1],[minB; 0; 0.5; rangepeak(1)],[maxB; 100*vmax; 3; rangepeak(2)],options);
peaks.B=x(1); peaks.A=x(2); peaks.sp=x(3); peaks.fp=x(4);
peaks.atot = peaks.B - peaks.C*log(peaks.fp) + peaks.A;


% 3 - compares to intertrial variations the rms noise
% [~,idpeak]=min(abs(peaks.fp-foi));
% rms=sqrt(median( (log(data(:,idpeak))-logpf(idpeak)).^2));
rms=(log(data)-repmat(logpf,Ntrials,1)).^2; rms=sqrt(median(rms(:)));


if or(peaks.A < rms, or( peaks.fp-peaks.sp<min(foi)-0.5,  peaks.fp+peaks.sp/3>max(foi)))
    peakfittedok=false;
else
    peakfittedok=true;
end

