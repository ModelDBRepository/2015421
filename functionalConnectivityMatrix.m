
% This code computes a functional connectivity matrix with a maximum lag of
% 20ms between spikes for the cross-correlation.

% These analyses were performed for O'Neill et al., Time-dependent
% homeostatic mechanisms underlie BDNF action on neural circuitry. Comms
% Bio, 2023.

% This function was written by Erin D. Anderson and can be
% accessed at https://www.seas.upenn.edu/~molneuro/

% Last Updated: 11/14/2023

function [FC] = functionalConnectivityMatrix(spikeindexes,spiketimes,...
    SimTimeInSeconds,nNeurons,maxlag)

if nargin == 4
    maxlag = 20; % [ms]
elseif nargin > 5
    warning('Too many inputs')
end

binsize = 1e-3;
edges = 0:binsize:SimTimeInSeconds;

%% Generate functional connectivity matrix

spiketimesNeuronal = cell(nNeurons,1);
spiketimesNeuronalDiscretized = zeros(nNeurons,length(edges)-1);
FC = ones(nNeurons)*0.5; % ones because we're not doing the diagonal
for ii = 1:nNeurons
    spiketimesNeuronal{ii} = spiketimes(spikeindexes == ii);
    spiketimesNeuronalDiscretized(ii,:) = histcounts(spiketimesNeuronal{ii},edges);
end

for ii = 1:nNeurons
    spiketimesNeuronal1 = spiketimesNeuronalDiscretized(ii,:);
    if numel(spiketimesNeuronal1) > 0 % if there's no connection, then don't bother correlating it with anything
        parfor jj = (ii+1):nNeurons
            spiketimesNeuronal2 = spiketimesNeuronalDiscretized(jj,:);
            FC(ii,jj) = max(xcorr(spiketimesNeuronal1,spiketimesNeuronal2,maxlag,'coeff'));
        end
    end
end

tril_FC = triu(FC)';
FC = tril_FC + triu(FC);
FC(isnan(FC)) = 0; % if there's no connection, set it equal to 0

