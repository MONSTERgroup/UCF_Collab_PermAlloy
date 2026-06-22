% Generalized code for analyzing various grain/sub grain/twin quanities
clear; close all;
currentPath = pwd;
pname = fullfile(currentPath); %Add folder name, if relevant
addpath(genpath(pname))

dataPath = fullfile(pname, 'Data');
if ~isfolder(dataPath)
    mkdir(dataPath);
end

imgPath = fullfile(pname, 'Images');
comImgPath = fullfile(imgPath, 'Combined Images');

if ~isfolder(imgPath)
    mkdir(imgPath);
    mkdir(comImgPath);
end

plotx2east
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','intoPlane'); % the BD would be Y
setMTEXpref('FontSize',14);

% Select desired group and initialize cell array
groupSets = struct();
groupSets.allXY = {...
    'As-Built 400mms 200W XY Scan 2_wSubGrains.mat',...
    'As-Built 900mms 200W XY Scan 2_wSubGrains.mat',...
    'As-Built 1200mms 200W XY Scan 2_wSubGrains.mat',...
    'HIP 120min 900mms 200W XY Scan 2_wSubGrains.mat',...
    'HIP 180min 900mms 200W XY Scan 2_wSubGrains.mat'};
groupSets.allXZ = {...
    'As-Built 400mms 200W XZ Scan 2_wSubGrains.mat',...
    'As-Built 900mms 200W XZ Scan 2_wSubGrains.mat',...
    'As-Built 1200mms 200W XZ Scan 4_wSubGrains.mat',...
    'HIP 120min 900mms 200W XZ Scan 2_wSubGrains.mat',...
    'HIP 180min 900mms 200W XZ Scan 2_wSubGrains.mat'};
groupSets.onlyHIP = {...
    'HIP 120min 900mms 200W XY Scan 2_wSubGrains.mat',...
    'HIP 180min 900mms 200W XY Scan 2_wSubGrains.mat',...
    'HIP 120min 900mms 200W XZ Scan 2_wSubGrains.mat',...
    'HIP 180min 900mms 200W XZ Scan 2_wSubGrains.mat'};
groupSets.onlyAsBuilt = {...
    'As-Built 400mms 200W XY Scan 2_wSubGrains.mat',...
    'As-Built 900mms 200W XY Scan 2_wSubGrains.mat',...
    'As-Built 1200mms 200W XY Scan 2_wSubGrains.mat',...    
    'As-Built 400mms 200W XZ Scan 2_wSubGrains.mat',...
    'As-Built 900mms 200W XZ Scan 2_wSubGrains.mat',...
    'As-Built 1200mms 200W XZ Scan 4_wSubGrains.mat'};
groupSets.compTopVsCenter = {...
    'As-Built 400mms 200W XZ Scan 2_wSubGrains.mat',....
    'As-Built 400mms 200W XZ top Scan 4_wSubGrains.mat',...
    'As-Built 900mms 200W XZ Scan 2_wSubGrains.mat',...
    'As-Built 900mms 200W XZ top Scan 1_wSubGrains.mat'};

chosenGroup = 'allXY';    % <-- change this string to pick a different group
groupToUse = groupSets.(chosenGroup);

%% Running twin analysis and loading in relevant file info for grain analysis
TwinStats = cell(length(groupToUse), 1);
allGrainData = cell(length(groupToUse), 1);

plotTwinFigures = false;     % <-- change to plot twin figures or not

searchPattern_GA = fullfile(dataPath, ['GrainAnalysis_' chosenGroup '.mat']);
searchPattern_TA = fullfile(dataPath, ['TwinAnalysis_' chosenGroup '.mat']);
GA_exists = exist(searchPattern_GA, 'file');
TA_exists = exist(searchPattern_TA, 'file');

if GA_exists && TA_exists
    response = input('Do you want to rerun your analysis? (y/n): ', 's');
    if ~strcmpi(response, 'y')
        fprintf('Loading existing analysis data...\n');
        load(searchPattern_GA, 'allGrainData');
        load(searchPattern_TA, 'TwinStats');
        return;  % Exit the script
    else
        fprintf('Rerunning analysis...\n');
    end
else
    fprintf('Running analysis for "%s"...\n', chosenGroup);
end

for i = 1:length(groupToUse)
    fprintf('Processing file %d of %d: %s\n', i, length(groupToUse), groupToUse{i});

    FileData = load(fullfile(pname, groupToUse{i}));
    ebsd = FileData.ebsd;
    grains = FileData.grains; 
    CS = ebsd.CS;

    allGrainData{i}.fileInfo = FileData.thisFile;
    allGrainData{i}.sampleName = FileData.thisFile.sampleName;
    allGrainData{i}.ebsd = ebsd;
    allGrainData{i}.grains = grains;
    allGrainData{i}.GrainResults = FileData.GrainResults;

    twinning = orientation.byAxisAngle(Miller(1,1,1,CS), 60*degree, CS, CS);

    % Twin Analysis (just in case we wanna do something with twin for HIPed sample later)
    gB = grains.boundary('Face Centered Cubic','Face Centered Cubic');
    isTwinning = angle(gB.misorientation,twinning) < 8.66*degree;
    twinBoundary = gB(isTwinning);
    
    if ~isempty(twinBoundary)
        [mergedGrains,parentId] = merge(grains,twinBoundary);
        % merge twins for GSD
        outerBoundaries_id = any(mergedGrains.boundary.grainId==0,2);
        grain_id = mergedGrains.boundary(outerBoundaries_id).grainId;
        grain_id(grain_id==0) = [];
        mergedGrains(grain_id) = [];
        gArea = grains.area;
        twinsPerGrain = zeros(mergedGrains.length,1);
        isTwin = true(grains.length,1);
        
        for j = 1:length(mergedGrains)
            % get child ids
            childId = find(parentId==j);
            for k = 1:length(childId)
                if grains(childId(k)).phaseId == 1
                    childId(k) = [];
                end
            end

            if isempty(childId)
                continue;
            end

            % cluster grains of similar orientations
            [fId,~] = calcCluster(grains(childId).meanOrientation,'maxAngle',...
                15*degree,'method','hierarchical','silent');
            % compute area of each cluster
            clusterArea = accumarray(fId,gArea(childId));
            % label the grains of largest cluster as original grain
            [~,fParent] = max(clusterArea);
            isTwin(childId(fId==fParent)) = false;
            % add number of grains to the array
            twinsPerGrain(j) = length(grains(childId(isTwin(childId))));
        end
        
        % compute the twinning statistics
        xMax = max(ebsd.prop.x); % width of ebsd scan (x)
        yMax = max(ebsd.prop.y); % height of ebsd scan (y)
        ebsd_Area = xMax*yMax; % area um^2
    
        TBDensity_ebsdArea = sum(twinBoundary.segLength)/ebsd_Area; % units: um^-1
        TBLengthPerGrain = accumarray(twinBoundary.grainId(:), repmat(twinBoundary.segLength, 2, 1), [length(grains) 1]);
        TBLengthPerGrainArea = accumarray(twinBoundary.grainId(:), repmat(twinBoundary.segLength, 2, 1), [length(gArea) 1]);

        isFCCGrain = grains.phaseId == 2;
        grainsFCC = grains(isFCCGrain);
        isTwin = isTwin(isFCCGrain);
        TBLengthPerGrain = TBLengthPerGrain(isFCCGrain);
        TBLengthPerGrainArea = TBLengthPerGrainArea(isFCCGrain);

        numTwins = length(grainsFCC(isTwin));
        [~,~,twinWidths] = grainsFCC(isTwin).fitEllipse;
        twinAreaFraction = sum(area(grainsFCC(isTwin)))/sum(area(grainsFCC)) * 100;
        twinnedFraction = (sum(twinsPerGrain>0)/length(mergedGrains)) * 100;
        
        if plotTwinFigures
            figure;
            plot(ebsd, ebsd.orientations)
            hold on
            plot(grains.boundary, 'lineWidth',1)
            hold on
            plot(mergedGrains.boundary,'lineWidth',1.5)
            export_fig(char(fullfile(imgPath, sprintf('%s %s.png', allGrainData{i}.sampleName, ' mergedGrains'))), '-m2'); 
            close
    
            figure;
            plot(ebsd, ebsd.orientations)
            hold on
            plot(grains.boundary,'lineWidth',1)
            hold on
            plot(twinBoundary,'lineWidth',0.8,'lineColor','white')
            export_fig(char(fullfile(imgPath, sprintf('%s %s.png', allGrainData{i}.sampleName, ' twinBoundaries'))), '-m2'); 
            close
        end
    else
        TBDensity_ebsdArea = 0;
        TBLengthPerGrain = [];
        TBLengthPerGrainArea = [];
        twinAreaFraction = 0;
        twinsPerGrain = zeros(grainsFCC.length,1);
        numTwins = 0;
        twinnedFraction = 0;
        twinWidths = [];
    end
    
    % Store results for this file
    TwinStats{i}.filename = groupToUse{i};
    TwinStats{i}.mergedGrains = mergedGrains;
    TwinStats{i}.TBDensity_ebsdArea = TBDensity_ebsdArea;
    TwinStats{i}.TBLengthPerGrain = TBLengthPerGrain;
    TwinStats{i}.TBLengthPerGrainArea = TBLengthPerGrainArea;
    TwinStats{i}.Merged_grains_w_twins = twinsPerGrain(twinsPerGrain > 0);
    TwinStats{i}.Twin_grain_area_frac = twinAreaFraction;
    TwinStats{i}.twin_grain_aspect_ratio = grainsFCC(isTwin).aspectRatio;
    TwinStats{i}.twin_grain_fraction = numel(grainsFCC(isTwin).id)/numel(grainsFCC.id);
    TwinStats{i}.numTwins = numTwins;
    TwinStats{i}.twinnedFraction = twinnedFraction;
    TwinStats{i}.twinWidths = twinWidths;
    
    fprintf('Completed analysis for %s\n', groupToUse{i});
end

save(char(fullfile(dataPath, strcat('TwinAnalysis_', chosenGroup))),"TwinStats","-v7.3");
save(char(fullfile(dataPath, strcat('GrainAnalysis_', chosenGroup))),"allGrainData","-v7.3");

%% Adapting colors to match with chosenGroup

numFiles = length(groupToUse);
defaultColors = lines(5);
if strcmp(chosenGroup, 'allXY') || strcmp(chosenGroup, 'allXZ')
    colors = defaultColors(1:numFiles, :);
elseif strcmp(chosenGroup, 'onlyAsBuilt')
    colorIndices = repmat([1, 2, 3], 1, ceil(numFiles/3));
    colors = defaultColors(colorIndices(1:numFiles), :);
elseif strcmp(chosenGroup, 'onlyHIP')
    colorIndices = repmat([4, 5], 1, ceil(numFiles/2));
    colors = defaultColors(colorIndices(1:numFiles), :);
else
    colorIndices = repmat(1:5, 1, ceil(numFiles/5));
    colors = defaultColors(colorIndices(1:numFiles), :);
end

%% Plotting grouped data

% Segmenting file names to get labels for plots
fullSampleLabels = cell(numFiles,1);
shortenLabels = cell(numFiles,1);

for labs = 1:numFiles
    nameIdx = allGrainData{labs}.sampleName;
    fileIdx = allGrainData{labs}.fileInfo;
    trimmedName = extractBefore(nameIdx, ' X');
    fullSampleLabels{labs} = trimmedName;

    if strcmp(fileIdx.sampleType, 'AsBuilt')
        parts = regexp(trimmedName, '(\S+)\s+(.*)', 'tokens', 'once');
        shortenLabels{labs} = parts{2};
    elseif strcmp(fileIdx.sampleType, 'HIP')
        parts = regexp(trimmedName, '(.*?min)\s+(.*)', 'tokens', 'once');
        shortenLabels{labs} = parts{1};
    else
        shortenLabels{labs} = trimmedName;
    end
end

% Equivalent Diameter - FCC Grains
allGrainEquDia = [];
dataSize_ones = [];

for f = 1:numFiles
    grainEquDia_f = allGrainData{f}.GrainResults.equDia_FCC;
    allGrainEquDia = [allGrainEquDia; grainEquDia_f(:)];
    numGrains_ones = ones(size(grainEquDia_f));
    dataSize_ones = [dataSize_ones;f*numGrains_ones];
end

figure('Position',[100 100 1000 800])
lognormalBoxplot(allGrainEquDia,dataSize_ones,'Colors',colors,'width', 0.45)
ylabel('Equivalent Diameter (\mum)')
ylim([0 200])
set(gcf,'Color','white')
xticklabels(shortenLabels)
PrettyPlotsSingle_BoxPlot(colors)
applyDashedLinesToBoxplot(chosenGroup, numFiles);
export_fig(char(fullfile(comImgPath, sprintf('%s %s.png', 'EquDia_grains', chosenGroup))), '-m2'); 
close

% Aspect Ratio - FCC Grains
allGrainAspRatio = [];
dataSize_ones = [];

for f = 1:numFiles
    grainAspRatio_f = allGrainData{f}.GrainResults.aspRatio_FCC;
    allGrainAspRatio = [allGrainAspRatio; grainAspRatio_f(:)];
    numGrains_ones = ones(size(grainAspRatio_f));
    dataSize_ones = [dataSize_ones;f*numGrains_ones];
end

figure('Position',[100 100 1000 800])
boxplot(allGrainAspRatio,dataSize_ones,'Colors',colors,'width', 0.45)
ylabel('Aspect Ratio')
ylim([0 20])
set(gcf,'Color','white')
xticklabels(shortenLabels)
PrettyPlotsSingle_BoxPlot(colors)
applyDashedLinesToBoxplot(chosenGroup, numFiles);
export_fig(char(fullfile(comImgPath, sprintf('%s %s.png', 'AspRatio_grains', chosenGroup))), '-m2'); 
close

% Sphericity - FCC Grains
allGrainSpher = [];
dataSize_ones = [];

for f = 1:numFiles
    grainSpher_f = allGrainData{f}.GrainResults.sphericity_FCC;
    allGrainSpher = [allGrainSpher; grainSpher_f(:)];
    numGrains_ones = ones(size(grainSpher_f));
    dataSize_ones = [dataSize_ones;f*numGrains_ones];
end

figure('Position',[100 100 1000 800])
boxplot(allGrainSpher,dataSize_ones,'Colors',colors,'width', 0.45)
ylabel('Sphercity Parameter')
ylim([0 0.5])
set(gcf,'Color','white')
xticklabels(shortenLabels)
PrettyPlotsSingle_BoxPlot(colors)
applyDashedLinesToBoxplot(chosenGroup, numFiles);
export_fig(char(fullfile(comImgPath, sprintf('%s %s.png', 'Sphericity_grains', chosenGroup))), '-m2'); 
close

%% Grain-based plots with mergedGrains instead
MergedGrainsData = cell(numFiles, 1);

for mg = 1:numFiles
    mergedGrains_FCC = TwinStats{mg}.mergedGrains('Face Centered Cubic');
    MergedGrainsData{mg}.mergedGrains_FCC = mergedGrains_FCC;
    MergedGrainsData{mg}.Mg_ED_FCC = 2*mergedGrains_FCC.equivalentRadius;
    MergedGrainsData{mg}.Mg_AR_FCC = mergedGrains_FCC.aspectRatio;
    MergedGrainsData{mg}.Mg_sphericity = mergedGrains_FCC.area ./ mergedGrains_FCC.perimeter('withInclusion') ./ mergedGrains_FCC.equivalentRadius;
end

save(char(fullfile(dataPath, strcat('MergedGrainsAnalysis_', chosenGroup))),"MergedGrainsData","-v7.3");

allGrainEquDia = [];
dataSize_ones = [];

for f = 1:numFiles
    grainEquDia_f = MergedGrainsData{f}.Mg_ED_FCC;
    allGrainEquDia = [allGrainEquDia; grainEquDia_f(:)];
    numGrains_ones = ones(size(grainEquDia_f));
    dataSize_ones = [dataSize_ones;f*numGrains_ones];
end

figure('Position',[100 100 1000 800])
lognormalBoxplot(allGrainEquDia,dataSize_ones,'Colors',colors,'width', 0.45)
ylabel('Equivalent Diameter (\mum)')
ylim([0 400])
set(gcf,'Color','white')
xticklabels(shortenLabels)
PrettyPlotsSingle_BoxPlot(colors)
applyDashedLinesToBoxplot(chosenGroup, numFiles);
export_fig(char(fullfile(comImgPath, sprintf('%s %s.png', 'EquDia_mergedGrains', chosenGroup))), '-m2'); 
close

% Aspect Ratio - FCC Merged Grains
allGrainAspRatio = [];
dataSize_ones = [];

for f = 1:numFiles
    grainAspRatio_f = MergedGrainsData{f}.Mg_AR_FCC;
    allGrainAspRatio = [allGrainAspRatio; grainAspRatio_f(:)];
    numGrains_ones = ones(size(grainAspRatio_f));
    dataSize_ones = [dataSize_ones;f*numGrains_ones];
end

figure('Position',[100 100 1000 800])
boxplot(allGrainAspRatio,dataSize_ones,'Colors',colors,'width', 0.45)
ylabel('Aspect Ratio')
ylim([0 20])
set(gcf,'Color','white')
xticklabels(shortenLabels)
PrettyPlotsSingle_BoxPlot(colors)
applyDashedLinesToBoxplot(chosenGroup, numFiles);
export_fig(char(fullfile(comImgPath, sprintf('%s %s.png', 'AspRatio_mergedGrains', chosenGroup))), '-m2'); 
close

% Sphericity - FCC Merged Grains
allGrainSpher = [];
dataSize_ones = [];

for f = 1:numFiles
    grainSpher_f = MergedGrainsData{f}.Mg_sphericity;
    allGrainSpher = [allGrainSpher; grainSpher_f(:)];
    numGrains_ones = ones(size(grainSpher_f));
    dataSize_ones = [dataSize_ones;f*numGrains_ones];
end

figure('Position',[100 100 1000 800])
boxplot(allGrainSpher,dataSize_ones,'Colors',colors,'width', 0.45)
ylabel('Sphercity Parameter')
ylim([0 0.5])
set(gcf,'Color','white')
xticklabels(shortenLabels)
PrettyPlotsSingle_BoxPlot(colors)
applyDashedLinesToBoxplot(chosenGroup, numFiles);
export_fig(char(fullfile(comImgPath, sprintf('%s %s.png', 'Sphericity_mergedGrains', chosenGroup))), '-m2'); 
close
