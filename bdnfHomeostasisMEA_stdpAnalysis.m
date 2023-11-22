
% This code compares the STDP-dependent exc-exc connection weights in
% different epochs and conditions.

% These analyses were performed for O'Neill et al., Time-dependent
% homeostatic mechanisms underlie BDNF action on neural circuitry. Comms
% Bio, 2023.

% This function was written by Erin D. Anderson and can be
% accessed at https://www.seas.upenn.edu/~molneuro/

% Last Updated: 11/18/2023

% function stdpAnalysis = bdnfHomeostasisMEA_stdpAnalysis(folderName, overwrite)

% Collect files to run
fileNames = dir(folderName);
for jj = length(fileNames):-1:1 % get rid of all the irrelevant files
    if fileNames(jj).isdir
        fileNames(jj) = [];
    elseif ~strcmp(fileNames(jj).name(end-3:end),'.mat')
        fileNames(jj) = [];
    end
end

interpLength = 10000;

% Initialize table
columnNames = {'FileName','CellDensity','PlateWidth','nNeurons','BDNF1','BDNF2','InjuryFraction1','InjuryFraction2',...
    'nEEConnections','STDP_pre','STDP_post','CDF_f','CDF_x','CDF_xCommon','CDF_fInterp'};
columnTypes = [{'string'}; repmat({'double'},[8,1]); repmat({'cell'},[length(columnNames)-1-8,1])];

STDPTableMEA = table('Size',[length(fileNames) length(columnNames)],...
    'VariableNames',columnNames,'VariableTypes',columnTypes);

% Run analysis
for jj = 1:length(fileNames)
    mat = matfile(fullfile(folderName,fileNames(jj).name));
    if ismember('nRegions', who('-file',fullfile(folderName,fileNames(jj).name)))
        if ~ismember('STDPDataMEA',who('-file',fullfile(folderName,fileNames(jj).name))) || overwrite
            STDPDataMEA = struct;
            
            % split into time segments
            statemon_dt = mat.w_STDP_timesInSeconds(1,2) - mat.w_STDP_timesInSeconds(1,1);
            
            STDPDataMEA.STDP_startstopIndex.pre = [mat.PreInjurySimTimeInSeconds/statemon_dt, ...
                (mat.PreInjurySimTimeInSeconds+mat.InjurySimTimeInSeconds)/statemon_dt-1]; % pre
            STDPDataMEA.STDP_startstopIndex.post = [(mat.PreInjurySimTimeInSeconds + mat.InjurySimTimeInSeconds*2)/statemon_dt, ...
                (mat.PreInjurySimTimeInSeconds+mat.InjurySimTimeInSeconds*4)/statemon_dt-1]; % analysis
            
            % remove injured neurons' connections - only need to look at
            % exc neurons bc w_STDP is for ee connections
            STDPDataMEA.eeIndexes = [mat.con_ee_i',mat.con_ee_j']; % get the indexes
            injuredIndexes = ismember(STDPDataMEA.eeIndexes(:,1), mat.exc_injuredNeurons1) ...
                | ismember(STDPDataMEA.eeIndexes(:,2), mat.exc_injuredNeurons1); % make sure neither origin or target neuron is injured
            STDPDataMEA.w_STDP_minusInjured = mat.w_STDP_statemon;
            STDPDataMEA.w_STDP_minusInjured(injuredIndexes,:) = []; % get rid of injured neurons' connections from weights
            STDPDataMEA.eeIndexes(injuredIndexes,:) = []; % get rid of injured neurons' indexes in the list too
            
            STDPDataMEA.STDP.pre = mean(STDPDataMEA.w_STDP_minusInjured(:,STDPDataMEA.STDP_startstopIndex.pre(1,1):...
                STDPDataMEA.STDP_startstopIndex.pre(1,2)),2);
            STDPDataMEA.STDP.post = mean(STDPDataMEA.w_STDP_minusInjured(:,STDPDataMEA.STDP_startstopIndex.post(1,1):...
                STDPDataMEA.STDP_startstopIndex.post(1,2)),2);
            
            [STDPDataMEA.CDF.f,STDPDataMEA.CDF.x] = ecdf(STDPDataMEA.STDP.post);
        else
            STDPDataMEA = mat.STDPDataMEA;
        end
        
        STDPTableMEA.FileName(jj) = fileNames(jj).name;
        STDPTableMEA.CellDensity(jj) = mat.CellDensity;
        STDPTableMEA.RecordingWidth(jj) = mat.RecordingWidth;
        STDPTableMEA.nNeurons(jj) = double(mat.N);
        STDPTableMEA.BDNF1(jj) = mat.BDNF(1,1);
        STDPTableMEA.BDNF2(jj) = mat.BDNF(1,2);
        STDPTableMEA.InjuryFraction1(jj) = mat.glutamate(1,1);
        STDPTableMEA.InjuryFraction2(jj) = mat.glutamate(1,2);
        
        clearvars mat
        STDPTableMEA.nEEConnections(jj) = size(STDPDataMEA.eeIndexes,1);
        STDPTableMEA.STDP_pre{jj} = STDPDataMEA.STDP.pre;
        STDPTableMEA.STDP_post{jj} = STDPDataMEA.STDP.post;
        STDPTableMEA.CDF_f{jj} = STDPDataMEA.CDF.f;
        STDPTableMEA.CDF_x{jj} = STDPDataMEA.CDF.x;
        
        STDPTableMEA.CDF_xCommon{jj} = linspace(0,16,interpLength); % interpolate because different range of w_ees for each simulation
        STDPTableMEA.CDF_fInterp{jj} = interp1(STDPTableMEA.CDF_x{jj}(2:end),STDPTableMEA.CDF_f{jj}(2:end),STDPTableMEA.CDF_xCommon{jj}');
        
%         save(fullfile(folderName,fileNames(jj).name),'STDPDataMEA','-append');
        
        clearvars STDPDataMEA
    else
        warndlg('Make sure to convert to MEA prior to running analysis!')
    end
    
end

% % save spikeTable for later
datestamp = clock;
datestamp = [num2str(datestamp(1)), sprintf('%02d',datestamp(2)),sprintf('%02d',datestamp(3))];
save(fullfile(folderName,['stdpTableMEA_',datestamp,'.mat']),'STDPTableMEA');

%% Split into different conditions
[combos,~,ind] = unique([STDPTableMEA.BDNF1,STDPTableMEA.BDNF2,...
    STDPTableMEA.InjuryFraction1,STDPTableMEA.InjuryFraction2],'rows');
RowsForEachCombo = accumarray(ind, find(ind), [], @(rows){rows});
TempLengths = [];

for i = 1:length(RowsForEachCombo)
    TempLengths = [TempLengths length(RowsForEachCombo{i})];
end

nCombos = size(combos,1);
cdf_x = STDPTableMEA.CDF_xCommon{1};
cdf_f = cell(1,length(RowsForEachCombo));
STDP_pre = NaN(interpLength*max(TempLengths),nCombos);
STDP_post = NaN(interpLength*max(TempLengths),nCombos);

for ii = 1:length(RowsForEachCombo)
    RowsForEachCombo{ii} = sort(RowsForEachCombo{ii}); % order matters because matching
    nRowsTotal = 0;
    for jj = 1:numel(RowsForEachCombo{ii})
        cdf_f{ii}(jj,:) = STDPTableMEA.CDF_fInterp{RowsForEachCombo{ii}(jj)};
        nRowsToAdd = length(STDPTableMEA.STDP_post{RowsForEachCombo{ii}(jj)});
        STDP_pre(nRowsTotal+1:nRowsTotal+nRowsToAdd,ii) = STDPTableMEA.STDP_pre{RowsForEachCombo{ii}(jj)};
        STDP_post(nRowsTotal+1:nRowsTotal+nRowsToAdd,ii) = STDPTableMEA.STDP_post{RowsForEachCombo{ii}(jj)};
        nRowsTotal = nRowsTotal+nRowsToAdd;
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

columnNames = {'Name','cdf_x_post','cdf_f_post','STDP_pre','STDP_post'};
columnTypes = {'string','double','cell','double','double'};

stdpAnalysis = table('Size',[nCombos length(columnNames)],...
    'VariableNames',columnNames,'VariableTypes',columnTypes);

stdpAnalysis.Name = comboNames;
stdpAnalysis.cdf_x_post = repmat(cdf_x,[nCombos,1]);
stdpAnalysis.cdf_f_post = cdf_f';
stdpAnalysis.w_ee_pre = STDP_pre'; % w_ee
stdpAnalysis.w_ee_post = STDP_post'; % w_ee
