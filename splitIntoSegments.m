
% This code separates spiketimes and spikeindexes into the different
% analysis "epochs."

% These analyses were performed for O'Neill et al., Time-dependent
% homeostatic mechanisms underlie BDNF action on neural circuitry. Comms
% Bio, 2023.

% This function was written by Erin D. Anderson and can be
% accessed at https://www.seas.upenn.edu/~molneuro/

% Last Updated: 11/14/2023

function [segmentedSpikeTimes,segmentedSpikeIndexes] = ...
    splitIntoSegments(spiketimes,spikeindexes,stabilizationTime,segmentLengths)

segmentedSpikeTimes = cell(length(segmentLengths),1);
segmentedSpikeIndexes = cell(length(segmentLengths),1);

segmentLengths = [0,segmentLengths];

for ii = 2:length(segmentLengths)
    tempIndexes = find(spiketimes >= sum(segmentLengths(1:ii-1)) + stabilizationTime &...
        spiketimes < sum(segmentLengths(1:ii)) + stabilizationTime);
    segmentedSpikeTimes{ii-1} = spiketimes(tempIndexes);
    if ~isempty(spikeindexes)
        segmentedSpikeIndexes{ii-1} = spikeindexes(tempIndexes);
    else
        segmentedSpikeIndexes{ii-1} = [];
    end
end

end