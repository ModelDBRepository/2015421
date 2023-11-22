
% This function loads in the simulation files and, if you want to overwrite
% the previous MEA conversion or this is the first time you're running the
% script, combines neuronal spike trains into a pseudo-MEA recording by
% combining neurons (exc or inh) that are physically proximal to each
% other. See Supplementary Figure 7 for more details. This function depends
% on neuronToMEA.m. 

% These analyses were performed for O'Neill et al., Time-dependent 
% homeostatic mechanisms underlie BDNF action on neural circuitry. Comms 
% Bio, 2023.

% This function was written by Erin D. Anderson and can be
% accessed at https://www.seas.upenn.edu/~molneuro/

% Last Updated: 11/14/2023 

function convertToMEA(folderName,overwrite)

nRowsColumns = 8; % hard coded to match MEA

% Load data
fileNames = dir(folderName);
for jj = length(fileNames):-1:1 % get rid of all the irrelevant files
    if fileNames(jj).isdir
        fileNames(jj) = [];
    elseif ~strcmp(fileNames(jj).name(end-3:end),'.mat')
        fileNames(jj) = [];
    end
end

for jj = 1:length(fileNames)
    disp(num2str(jj))
    mat = matfile(fullfile(folderName,fileNames(jj).name));
    if ~ismember('nRegions',who('-file',fullfile(folderName,...
            fileNames(jj).name))) || overwrite
        
        % group together excitatory and inhibitory neurons
        spiketimes = double([mat.exc_spiketimes,mat.inh_spiketimes]);
        spikeindexes = double([mat.exc_spikeindexes,mat.inh_spikeindexes+mat.NE]);
        
        % convert to MEA
        [nRegions,spikeIndexesGrid,spikeTimesGrid,neuronsGrid] = ...
            neuronToMEA(mat.NeuronXPosition,mat.NeuronYPosition,...
            mat.RecordingWidth,spiketimes,spikeindexes,nRowsColumns);
        
        clear spiketimes spikeindexes
        
        % get rid of non-recording electrodes (hardcoded to match
        % actual MEA)
        electrodesToRemove = [1,8,25,56,64];
        
        % hard coded out of laziness (sorry)
        indexesOfElectrodesToRemove = spikeIndexesGrid == 1 | ...
            spikeIndexesGrid == 8 | spikeIndexesGrid == 25 | ...
            spikeIndexesGrid == 56 | spikeIndexesGrid == 64;
        
        spikeTimesGrid(indexesOfElectrodesToRemove) = [];
        spikeIndexesGrid(indexesOfElectrodesToRemove) = [];
        
        for kk = length(electrodesToRemove):-1:1
            spikeIndexesGrid(spikeIndexesGrid > electrodesToRemove(kk)) = ...
                spikeIndexesGrid(spikeIndexesGrid > electrodesToRemove(kk)) - 1;
        end
        
        nRegions = nRegions - length(electrodesToRemove);
        
%         save(fullfile(fileNames(jj).folder,fileNames(jj).name),'-v7.3') % if it's not already v7.3
        save(fullfile(fileNames(jj).folder,fileNames(jj).name),'nRegions','spikeIndexesGrid','spikeTimesGrid','neuronsGrid','-append');

    end
end
end
