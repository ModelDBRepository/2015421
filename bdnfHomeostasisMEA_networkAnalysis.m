
% This code computes the Burstlet Rate & Fano factor

% These analyses were performed for O'Neill et al., Time-dependent
% homeostatic mechanisms underlie BDNF action on neural circuitry. Comms
% Bio, 2023.

% This function was written by Erin D. Anderson and can be
% accessed at https://www.seas.upenn.edu/~molneuro/

% Last Updated: 11/14/2023

function networkAnalysis = bdnfHomeostasisMEA_networkAnalysis(folderName, overwrite)

% Collect files to run
fileNames = dir(folderName);
for jj = length(fileNames):-1:1 % get rid of all the irrelevant files
    if fileNames(jj).isdir
        fileNames(jj) = [];
    elseif ~strcmp(fileNames(jj).name(end-3:end),'.mat')
        fileNames(jj) = [];
    end
end

% Initialize table
columnNames = {'FileName','CellDensity','RecordingWidth','nNeurons','nRegions','BDNF1','BDNF2','InjuryFraction1','InjuryFraction2',...
    'FC_Pre','FC_Post',...
    'LocalEfficiency_Pre','LocalEfficiency_Post','LocalEfficiency_Norm','LocalEfficiency_Delta'};
columnTypes = [{'string'}; repmat({'double'},[8,1]);
    repmat({'cell'},[2,1]); repmat({'cell'},[4,1])];

networkTableMEA = table('Size',[length(fileNames) length(columnNames)],...
    'VariableNames',columnNames,'VariableTypes',columnTypes);

% Run analysis
for jj = 1:length(fileNames)
    mat = matfile(fullfile(folderName,fileNames(jj).name));
    if ismember('nRegions', who('-file',fullfile(folderName,fileNames(jj).name)))
        if ~ismember('networkDataMEA',who('-file',fullfile(folderName,...
                fileNames(jj).name))) || overwrite
            networkDataMEA = struct;
            
            % split into pre, during, and post segments
            [segmentedSpikeTimes,segmentedSpikeIndexes] = ...
                splitIntoSegments(mat.spikeTimesGrid,mat.spikeIndexesGrid,...
                mat.PreInjurySimTimeInSeconds,...
                [mat.SimTimeInSeconds,mat.InjurySimTimeInSeconds,mat.InjurySimTimeInSeconds,mat.PostInjurySimTimeInSeconds]);
            
            % FC
            networkDataMEA.FC.pre = functionalConnectivityMatrix(segmentedSpikeIndexes{1},...
                segmentedSpikeTimes{1}-mat.PreInjurySimTimeInSeconds,mat.SimTimeInSeconds,mat.nRegions);
            networkDataMEA.FC.post = functionalConnectivityMatrix(segmentedSpikeIndexes{4},...
                segmentedSpikeTimes{4}-mat.PreInjurySimTimeInSeconds-mat.SimTimeInSeconds-mat.InjurySimTimeInSeconds*2,...
                mat.PostInjurySimTimeInSeconds,mat.nRegions);
            
            % LE
            networkDataMEA.LE.pre = efficiency_wei(networkDataMEA.FC.pre,2);
            networkDataMEA.LE.post = efficiency_wei(networkDataMEA.FC.post,2);
            networkDataMEA.LE.norm = networkDataMEA.LE.post./networkDataMEA.LE.pre;
            networkDataMEA.LE.delta = networkDataMEA.LE.post - networkDataMEA.LE.pre;
            
%             save(fullfile(fileNames(jj).folder,fileNames(jj).name),'networkDataMEA','-append');
        else
            networkDataMEA = mat.networkDataMEA;
        end
               
        networkTableMEA.FileName(jj) = fileNames(jj).name;
        networkTableMEA.CellDensity(jj) = mat.CellDensity;
        networkTableMEA.RecordingWidth(jj) = mat.RecordingWidth;
        networkTableMEA.nNeurons(jj) = double(mat.N);
        networkTableMEA.nRegions(jj) = mat.nRegions;
        networkTableMEA.BDNF1(jj) = mat.BDNF(1,1);
        networkTableMEA.BDNF2(jj) = mat.BDNF(1,2);
        networkTableMEA.InjuryFraction1(jj) = mat.glutamate(1,1);
        networkTableMEA.InjuryFraction2(jj) = mat.glutamate(1,2);
        
        networkTableMEA.FC_Pre{jj} = networkDataMEA.FC.pre;
        networkTableMEA.FC_Post{jj} = networkDataMEA.FC.post;
        
        networkTableMEA.LocalEfficiency_Pre{jj} = networkDataMEA.LE.pre;
        networkTableMEA.LocalEfficiency_Post{jj} = networkDataMEA.LE.post;
        networkTableMEA.LocalEfficiency_Norm{jj} = networkDataMEA.LE.norm;
        networkTableMEA.LocalEfficiency_Delta{jj} = networkDataMEA.LE.delta;
        
        clearvars networkDataMEA
    else
        warndlg('Make sure to convert to MEA prior to running analysis!')
    end
end

% % save networkTable for later
datestamp = clock;
datestamp = [num2str(datestamp(1)), sprintf('%02d',datestamp(2)),sprintf('%02d',datestamp(3))];
save(fullfile(folderName,['networkTableMEA_',datestamp,'.mat']),'networkTableMEA');

%% Split into different conditions
[combos,~,jj] = unique([networkTableMEA.BDNF1,networkTableMEA.BDNF2,...
    networkTableMEA.InjuryFraction1,networkTableMEA.InjuryFraction2],'rows');
RowsForEachCombo = accumarray(jj, find(jj), [], @(rows){rows});
TempLengths = [];

for i = 1:length(RowsForEachCombo)
    TempLengths = [TempLengths length(RowsForEachCombo{i})];
end

nCombos = size(combos,1);
nRegions = networkTableMEA.nRegions(1);

LEpre = NaN(max(TempLengths)*nRegions,nCombos);
LEpost = NaN(max(TempLengths)*nRegions,nCombos);

for ii = 1:length(RowsForEachCombo)
    RowsForEachCombo{ii} = sort(RowsForEachCombo{ii}); % order matters because matching
    for jj = 1:numel(RowsForEachCombo{ii})
        LEpre(nRegions*(jj-1)+1:nRegions*jj,ii) = networkTableMEA.LocalEfficiency_Pre{RowsForEachCombo{ii}(jj)};
        LEpost(nRegions*(jj-1)+1:nRegions*jj,ii) = networkTableMEA.LocalEfficiency_Post{RowsForEachCombo{ii}(jj)};
    end
end

%% Make Summary Table

comboNames = cell(nCombos,1);
for ii = 1:nCombos
    if combos(ii,1) == 0 && combos(ii,2) == 0 && ... % no BDNF
        combos(ii,3) == 0 && combos(ii,4) == 0 % or injury
            comboNames{ii} = 'Control';
    elseif combos(ii,3) ~= 0 % if injured
        if combos(ii,2) ~= 0 % and BDNF'd
            comboNames{ii} = 'Injury then BDNF';
        else % no BDNF
            comboNames{ii} = 'Injury';
        end
    end
end

columnNames = {'Name','LocalEfficiencyPre','LocalEfficiencyPost'};
columnTypes = {'string','double','double'};

networkAnalysis = table('Size',[nCombos length(columnNames)],...
    'VariableNames',columnNames,'VariableTypes',columnTypes);

networkAnalysis.Name = comboNames;
networkAnalysis.LocalEfficiencyPre = LEpre';
networkAnalysis.LocalEfficiencyPost = LEpost';
