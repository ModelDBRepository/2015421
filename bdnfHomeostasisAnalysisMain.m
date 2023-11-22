
% This script calls the functions required to fully analyze the data output
% from Python - bdnfHomeostasisSimulation.py. See the individual functions 
% called here for more detailed information. These analyses were performed 
% for O'Neill et al., Time-dependent homeostatic mechanisms underlie BDNF 
% action on neural circuitry. Comms Bio, 2023.

% The functions called here were written by Erin D. Anderson and can be 
% accessed at https://www.seas.upenn.edu/~molneuro/

% Last Updated: 11/18/2023

%% Analysis
folderName = uigetdir; % folder where BDNF simulation data is stored
overwrite = false; % whether to overwrite data structs in simulation mat 
                   % files, e.g. if the analysis functions changed

convertToMEA(folderName, overwrite); % converts individual neuron recording to MEA recording
                                     % writes directly into file
                                     
spikeAnalysis = bdnfHomeostasisMEA_spikeAnalysis(folderName, overwrite); 
spikeAnalysisStruct = table2struct(spikeAnalysis);

networkAnalysis = bdnfHomeostasisMEA_networkAnalysis(folderName, overwrite);
networkAnalysisStruct = table2struct(networkAnalysis);

stdpAnalysis = bdnfHomeostasisMEA_stdpAnalysis(folderName, overwrite);
stdpAnalysisStruct = table2struct(stdpAnalysis);



