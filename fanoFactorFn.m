
% This function calculates the Fano factor based on a 100ms windowsize. 

% This analysis was performed for O'Neill et al., Time-dependent 
% homeostatic mechanisms underlie BDNF action on neural circuitry. Comms 
% Bio, 2023.

% This function was written by Erin D. Anderson and can be
% accessed at https://www.seas.upenn.edu/~molneuro/

% Last Updated: 11/14/2023 

function FF = fanoFactorFn(spikeindexes,spiketimes,SimTimeInSeconds,PreInjurySimTime,nNeurons)

windowSize = 0.1; % [s] 100 ms
nWindows = SimTimeInSeconds/windowSize; % divide into 100 ms intervals
FF = NaN(nNeurons,1);
nSpikes = NaN(nWindows,nNeurons);

for ii = 1:nNeurons
    neuronalSpikes = spiketimes(spikeindexes == ii-1); % get the spikes that correspond to current neuron (ii-1 to account for differential indexing in Python vs Matlab)
    for jj = 1:nWindows
        nSpikes(ii,jj) = nansum(neuronalSpikes > ((jj-1)*windowSize + PreInjurySimTime) & ...
            neuronalSpikes <= ((jj)*windowSize + PreInjurySimTime)); % count spikes in that interval after PreInjurySimTime
    end
    FF(ii) = nanvar(nSpikes(ii,:))/nanmean(nSpikes(ii,:)); % Fano Factor is the ratio of the sample
            % variance to the sample mean
end

end