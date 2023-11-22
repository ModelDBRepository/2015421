
% This function does the actual work of combining the neurons in a
% spatially-dependent manner to simulate an MEA recording. 

% These analyses were performed for O'Neill et al., Time-dependent 
% homeostatic mechanisms underlie BDNF action on neural circuitry. Comms 
% Bio, 2023.

% This function was written by Erin D. Anderson and can be
% accessed at https://www.seas.upenn.edu/~molneuro/

% Last Updated: 11/14/2023 

function [nRegions,spikeIndexesGrid,spikeTimesGrid,neuronsGrid] = ...
    neuronToMEA(NeuronXPosition,NeuronYPosition,RecordingWidth,spiketimes,spikeindexes,nRowsColumns)

increments = linspace(0,RecordingWidth,nRowsColumns+1); % divide plate into rows and columns

nRegions = nRowsColumns^2;
spikeIndexesGrid = []; % spike indexes based on regions in the grid
spikeTimesGrid = []; % spike times based on regions in the grid
neuronsGrid = cell(nRegions,1); % identity of the neurons in each region in the grid
spikeTimesThisRegion = cell(nRegions,1);
for xx = 1:nRowsColumns
    neuronsY = NeuronYPosition > increments(xx) & NeuronYPosition <= increments(xx+1); % neurons in correct y range
    for yy = 1:nRowsColumns
        neuronsX = NeuronXPosition > increments(yy) & NeuronXPosition <= increments(yy+1); % neurons in correct x range
        thisRegion = nRowsColumns*(xx-1)+yy; % convert 2-d region to 1-d
        neuronsGrid{thisRegion} = find(neuronsX & neuronsY) - 1; % subtract 1 to deal with Python indexing
        
        % split out the spike times for the neurons in each region
        spikeTimesThisRegion{thisRegion} = [];
        for ii = 1:numel(neuronsGrid{thisRegion})
            spikeTimesThisRegion{thisRegion} = [spikeTimesThisRegion{thisRegion}, ...
                spiketimes(spikeindexes == neuronsGrid{thisRegion}(ii))];
        end
              
        % convert spiketimes, spikeindexes to grid region activity instead
        % of neuron activity
        spikeTimesGrid = [spikeTimesGrid spikeTimesThisRegion{thisRegion}];
        spikeIndexesGrid = [spikeIndexesGrid repmat(thisRegion,[1,length(spikeTimesThisRegion{thisRegion})])];
    end
end

end
