
% This function was originally written for Rodriguez, et al., Cytosolic PSD-95 
% Interactor Regulates AMPA Receptor Signaling Independent of PSD-95 
% Binding. Network Neurosci, 2020.

% This function takes in the spiketimes and the length of the simulation to
% calculate the "network spike" features: the average inter-network-spike
% interval, the standard dev of the INSI, the number of network spikes for
% the simulation, the network spikes per second, and the coefficient of
% variation for the network spikes.

% This function was written by Erin D. Anderson and can be
% accessed at https://www.seas.upenn.edu/~molneuro/

% Last Updated: 11/14/2023

function [insi_avg,insi_std,NScount,NSrate,CV,nSpikesPerBurst_avg,nSpikesPerBurst_std,burstWidthAvg,burstWidthStd] = ...
    INSIbyElectrode(spikeindexes,spiketimes,SimTimeInSeconds,StabilizationTime,nUnits)

if nargin == 2
    StabilizationTime = 10;
    warning('Setting Stabilization Time to 10 s')
end

thr = 4; % need at least 4 events/bin to call it a network spike
timeBin = 100e-3; % [s] - thr and timeBin derived from Rodriguez thesis

insi_avg = zeros(nUnits,1);
insi_std = zeros(nUnits,1);
NScount = zeros(nUnits,1);
NSrate = zeros(nUnits,1);
CV = zeros(nUnits,1);
nSpikesPerBurst_avg = zeros(nUnits,1);
nSpikesPerBurst_std = zeros(nUnits,1);
burstWidthAvg = zeros(nUnits,1);
burstWidthStd = zeros(nUnits,1);

for ii = 1:nUnits    
    
    spiketimes_post = spiketimes(spiketimes(spikeindexes == ii) > StabilizationTime); % spiketimes that happen
    % after the settling period for this particular neuron
    
    if numel(spiketimes_post) > 1 && spiketimes_post(end)-spiketimes_post(1) > timeBin
        % divide the spiketimes into bins to see if reach threshold for network spike
        h = histcounts(spiketimes_post,'BinEdges',...
            [spiketimes_post(1):timeBin:spiketimes_post(end)]);
        rate = h/60.0/timeBin; % h[0] is the counts and h[1] is the bin edges
        i_start = find((rate(2:end)>=thr).*(rate(1:end-1)<thr)>0); % these are
                                                                    % indexes in the
                                                                    % histcounts
                                                                    % bins
                       % find all the indices > threshold, then check the index
                       % before, make sure it's 0, then multiply to get the start
                       % index
        i_end = find((rate(1:end-1)>=thr).*(rate(2:end)<thr)>0);
        
        if length(i_end) > length(i_start)
            i_end = i_end(2:end); % in case it thinks index = 1 is the end of a burst
        elseif length(i_start) > length(i_end)
            i_end(end+1) = i_start(end); % in case it thinks the last index is the start of a burst
        end
        
        if numel(i_start) <= 1
            insi = NaN; % if there are no network spikes
            burstWidth = NaN;
        else
            insi = ( i_start(2:end) - i_start(1:end-1) )*timeBin; % convert to time intervals
            burstWidth = (i_end-i_start)*timeBin;
        end
        insi_avg(ii) = mean(insi);
        insi_std(ii) = std(insi);
        NScount(ii) = length(i_start);
        NSrate(ii) = NScount(ii)/SimTimeInSeconds;
        CV(ii) = insi_std(ii)/insi_avg(ii);
        burstWidthAvg(ii) = mean(burstWidth);
        burstWidthStd(ii) = std(burstWidth);
        
        if NScount(ii) ~=0
            nSpikesPerBurst_avg(ii) = mean(h(i_start+1)); % i_start is the onset and i_start+1 is peak
            nSpikesPerBurst_std(ii) = std(h(i_start+1));
        else
            nSpikesPerBurst_avg(ii) = 0;
            nSpikesPerBurst_std(ii) = 0;
        end
    else
        insi_avg(ii) = 0;
        insi_std(ii) = 0;
        NScount(ii) = 0;
        NSrate(ii) = 0;
        CV(ii) = 0;
        nSpikesPerBurst_avg(ii) = 0;
        nSpikesPerBurst_std(ii) = 0;
        burstWidthAvg(ii) = 0;
        burstWidthStd(ii) = 0;
    end
end

end