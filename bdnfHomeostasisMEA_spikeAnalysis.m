
% This code computes the Burstlet Rate & Fano factor 

% These analyses were performed for O'Neill et al., Time-dependent
% homeostatic mechanisms underlie BDNF action on neural circuitry. Comms
% Bio, 2023.

% This function was written by Erin D. Anderson and can be
% accessed at https://www.seas.upenn.edu/~molneuro/

% Last Updated: 11/14/2023

function spikeAnalysis = bdnfHomeostasisMEA_spikeAnalysis(folderName, overwrite)

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
columnNames = {'FileName','CellDensity','RecordingWidth','nNeurons','nRegions',...
    'BDNF1','BDNF2','InjuryFraction1','InjuryFraction2'...
    'BurstRate_Pre','BurstRate_During','BurstRate_Post','BurstRate_Norm','BurstRate_Delta',...
    'FF_Pre','FF_During1','FF_During2','FF_Post','FF_Norm','FF_Delta'};
columnTypes = [{'string'}; repmat({'double'},[8,1]); repmat({'cell'},[length(columnNames)-1-8,1])];

spikeTableMEA = table('Size',[length(fileNames) length(columnNames)],...
    'VariableNames',columnNames,'VariableTypes',columnTypes);

% Run analysis
for jj = 1:length(fileNames)
    mat = matfile(fullfile(folderName,fileNames(jj).name));
    if ismember('nRegions', who('-file',fullfile(folderName,fileNames(jj).name)))
        if ~ismember('spikeDataMEA',who('-file',fullfile(folderName,fileNames(jj).name))) || overwrite
            spikeDataMEA = struct;
            
            % Split into pre, during, and post segments
            [segmentedSpikeTimes,segmentedSpikeIndexes] = ...
                splitIntoSegments(mat.spikeTimesGrid,mat.spikeIndexesGrid,...
                mat.PreInjurySimTimeInSeconds,...
                [mat.SimTimeInSeconds,mat.InjurySimTimeInSeconds,mat.InjurySimTimeInSeconds,mat.PostInjurySimTimeInSeconds]);
            
            % Burst rate
            [~,~,~,spikeDataMEA.burstRate.pre,~,~,...
                spikeDataMEA.nSpikesPerBurst_std.pre] = INSIbyElectrode(segmentedSpikeIndexes{1},segmentedSpikeTimes{1},...
                mat.SimTimeInSeconds,0,mat.nRegions);
            [~,~,~,spikeDataMEA.burstRate.during1,~,~,...
                spikeDataMEA.nSpikesPerBurst_std.during1] = INSIbyElectrode(segmentedSpikeIndexes{2},segmentedSpikeTimes{2},...
                mat.InjurySimTimeInSeconds,0,mat.nRegions);
            [~,~,~,spikeDataMEA.burstRate.during2,~,~,....
                spikeDataMEA.nSpikesPerBurst_std.during2] = INSIbyElectrode(segmentedSpikeIndexes{3},segmentedSpikeTimes{3},...
                mat.InjurySimTimeInSeconds,0,mat.nRegions);
            [~,~,~,spikeDataMEA.burstRate.post,~,~,....
                spikeDataMEA.nSpikesPerBurst_std.post] = INSIbyElectrode(segmentedSpikeIndexes{4},segmentedSpikeTimes{4},...
                mat.PostInjurySimTimeInSeconds,0,mat.nRegions);
            spikeDataMEA.burstRate.norm = spikeDataMEA.burstRate.post./spikeDataMEA.burstRate.pre;
            spikeDataMEA.burstRate.delta = spikeDataMEA.burstRate.post - spikeDataMEA.burstRate.pre;
            
            
            % Fano factor
            spikeDataMEA.FF.pre = fanoFactorFn(segmentedSpikeIndexes{1}-1,segmentedSpikeTimes{1},mat.SimTimeInSeconds,...
                mat.PreInjurySimTimeInSeconds,mat.nRegions); % have to manually input pre-sim time for FF bc need to know exact start time
            spikeDataMEA.FF.during1 = fanoFactorFn(segmentedSpikeIndexes{2}-1,segmentedSpikeTimes{2},mat.InjurySimTimeInSeconds,...
                mat.PreInjurySimTimeInSeconds + mat.SimTimeInSeconds,mat.nRegions);
            spikeDataMEA.FF.during2 = fanoFactorFn(segmentedSpikeIndexes{3}-1,segmentedSpikeTimes{3},mat.InjurySimTimeInSeconds,...
                mat.PreInjurySimTimeInSeconds + mat.SimTimeInSeconds + mat.InjurySimTimeInSeconds,mat.nRegions);
            spikeDataMEA.FF.post = fanoFactorFn(segmentedSpikeIndexes{4}-1,segmentedSpikeTimes{4},mat.PostInjurySimTimeInSeconds,...
                mat.PreInjurySimTimeInSeconds + mat.SimTimeInSeconds + mat.InjurySimTimeInSeconds*2,mat.nRegions);
            spikeDataMEA.FF.norm = spikeDataMEA.FF.post./spikeDataMEA.FF.pre;
            spikeDataMEA.FF.delta = spikeDataMEA.FF.post - spikeDataMEA.FF.pre;
            
%             save(fullfile(folderName,fileNames(jj).name),'spikeDataMEA','-append');
        else
            spikeDataMEA = mat.spikeDataMEA;
        end
        
        spikeTableMEA.FileName(jj) = fileNames(jj).name;
        spikeTableMEA.CellDensity(jj) = mat.CellDensity;
        spikeTableMEA.RecordingWidth(jj) = mat.RecordingWidth;
        spikeTableMEA.nNeurons(jj) = double(mat.N);
        spikeTableMEA.nRegions(jj) = mat.nRegions;
        spikeTableMEA.BDNF1(jj) = mat.BDNF(1,1);
        spikeTableMEA.BDNF2(jj) = mat.BDNF(1,2);
        spikeTableMEA.InjuryFraction1(jj) = mat.glutamate(1,1);
        spikeTableMEA.InjuryFraction2(jj) = mat.glutamate(1,2);
        
        spikeTableMEA.BurstRate_Pre{jj} = spikeDataMEA.burstRate.pre;
        spikeTableMEA.BurstRate_During1{jj} = spikeDataMEA.burstRate.during1;
        spikeTableMEA.BurstRate_During2{jj} = spikeDataMEA.burstRate.during2;
        spikeTableMEA.BurstRate_Post{jj} = spikeDataMEA.burstRate.post;
        spikeTableMEA.BurstRate_Norm{jj} = spikeDataMEA.burstRate.norm;
        spikeTableMEA.BurstRate_Delta{jj} = spikeDataMEA.burstRate.delta;
        
        spikeTableMEA.FF_Pre{jj} = spikeDataMEA.FF.pre;
        spikeTableMEA.FF_During1{jj} = spikeDataMEA.FF.during1;
        spikeTableMEA.FF_During2{jj} = spikeDataMEA.FF.during2;
        spikeTableMEA.FF_Post{jj} = spikeDataMEA.FF.post;
        spikeTableMEA.FF_Norm{jj} = spikeDataMEA.FF.norm;
        spikeTableMEA.FF_Delta{jj} = spikeDataMEA.FF.delta;
        
        clearvars spikeDataMEA
    else
        warndlg('Make sure to convert to MEA prior to running analysis!')
    end
end

% % save spikeTable for later
datestamp = clock;
datestamp = [num2str(datestamp(1)), sprintf('%02d',datestamp(2)),sprintf('%02d',datestamp(3))];
save(fullfile(folderName,['spikeTableMEA_',datestamp,'.mat']),'spikeTableMEA');

%% Split into different conditions
[combos,~,jj] = unique([spikeTableMEA.BDNF1,spikeTableMEA.BDNF2,...
    spikeTableMEA.InjuryFraction1,spikeTableMEA.InjuryFraction2],'rows');
RowsForEachCombo = accumarray(jj, find(jj), [], @(rows){rows});
TempLengths = [];

for i = 1:length(RowsForEachCombo)
    TempLengths = [TempLengths length(RowsForEachCombo{i})];
end

nCombos = size(combos,1);
nRegions = spikeTableMEA.nRegions(1); % they all have the same number of regions

burstRatePre = NaN(max(TempLengths)*nRegions,nCombos);
burstRatePost = NaN(max(TempLengths)*nRegions,nCombos);
fanoFactorPre = NaN(max(TempLengths)*nRegions,nCombos);
fanoFactorPost = NaN(max(TempLengths)*nRegions,nCombos);

for ii = 1:length(RowsForEachCombo)
    RowsForEachCombo{ii} = sort(RowsForEachCombo{ii}); % order matters because matching
    for jj = 1:numel(RowsForEachCombo{ii})
        burstRatePre(nRegions*(jj-1)+1:nRegions*jj,ii) = spikeTableMEA.BurstRate_Pre{RowsForEachCombo{ii}(jj)};
        burstRatePost(nRegions*(jj-1)+1:nRegions*jj,ii) = spikeTableMEA.BurstRate_Post{RowsForEachCombo{ii}(jj)};
        fanoFactorPre(nRegions*(jj-1)+1:nRegions*jj,ii) = spikeTableMEA.FF_Pre{RowsForEachCombo{ii}(jj)};
        fanoFactorPost(nRegions*(jj-1)+1:nRegions*jj,ii) = spikeTableMEA.FF_Post{RowsForEachCombo{ii}(jj)};
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

columnNames = {'Name','BurstRatePre','BurstRatePost','FanoFactorPre','FanoFactorPost'};
columnTypes = {'string','double','double','double','double'};

spikeAnalysis = table('Size',[nCombos length(columnNames)],...
    'VariableNames',columnNames,'VariableTypes',columnTypes);

spikeAnalysis.Name = comboNames;
spikeAnalysis.BurstRatePre = burstRatePre';
spikeAnalysis.BurstRatePost = burstRatePost';
spikeAnalysis.FanoFactorPre = fanoFactorPre';
spikeAnalysis.FanoFactorPost = fanoFactorPost';


