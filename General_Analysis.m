%% UCF Collab EBSD Analysis

clear; close all

% Setting up MTEX related info for plotting
plotx2east
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','intoPlane');
setMTEXpref('FontSize',14);

% crystal symmetry
CS = {... 
  'notIndexed',...
  crystalSymmetry('432', [4 4 4], 'mineral', 'Face Centered Cubic', 'color', [0.53 0.81 0.98])};
csFCC = CS{2}.Laue;

% extracting all relevant files
pname = pwd; % should be in folder with all the .ang files
addpath(genpath(pname))

angFiles = dir(fullfile(pname,'*.ang'));
allFiles = {angFiles.name}';

matFiles = dir(fullfile(pname,'*.mat'));
matNames = {matFiles.name}';
allMatBases = cellfun(@(n) erase(n,'.mat'),matNames,'UniformOutput',false);

% storing relevant info for later processing/analysis
allFileInfo = {};

for i = 1:length(allFiles)
    trimmed  = regexprep(allFiles{i}, '^\d{2}[\.-]\d{2}[\.-]\d{2}\s*', '');
    samp_name = regexprep(trimmed, '^[^\s]+\s+', '');
    samp_name = erase(samp_name, '.ang');
    firstPortion = extractBefore(samp_name, ' ');

    match = regexp(samp_name, ...
        '(?<params>\d+mms\s+\d+W)\s+(?<dir>[A-Z]{2})\s+Scan\s+(?<scan>\d+)', ...
        'names');

    thisFile = struct( ...
        'filename',        allFiles{i}, ...
        'sampleName',      samp_name, ...
        'sampleType',      firstPortion, ...
        'buildParams',     match.params, ...
        'buildDirection',  match.dir, ...
        'scanNumber',      str2double(match.scan));

    if strcmpi(firstPortion, 'HIP')
        htInfo = regexp(samp_name, '(\d+)\s*min', 'tokens', 'once');
        thisFile.htDurationMin = str2double(htInfo{1});
    else
        thisFile.htDurationMin = [];   % empty for As-Built
    end

    allFileInfo{end+1} = thisFile;
end

% initial ebsd processing (only need to do once for each .ang file)
for j = 1:numel(allFileInfo)
    samp_name = allFileInfo{j}.sampleName;
    
    if any(strcmpi(samp_name, allMatBases))
        fprintf('Data has been processed already for %s. Skipping...\n', allFiles{j});
    else
        fname = [pname filesep allFiles{j}];
        ebsd = EBSD.load(fname,csFCC,'interface','ang', 'convertEuler2SpatialReferenceFrame','setting 2'); 
        if strcmpi(allFileInfo{j}.buildDirection, 'XZ')
            ebsd = rotate(ebsd,180*degree,'keepEuler'); % rotates build direction to align with y-axis
        end
        ebsd_raw = ebsd; % making copy, just in case

        % plotting pre-processed EBSD IPF and IQ maps
        figure;
        plot(ebsd,ebsd.orientations) %,'micronbar','off') Uncomment, if desired
        export_fig(sprintf('%s %s.png', samp_name, ' rawEBSD_ipf'), '-m2');
        close

        figure;
        plot(ebsd,ebsd.iq) %,'micronbar','off') Uncomment, if desired
        colormap('gray')
        export_fig(sprintf('%s %s.png', samp_name, ' rawEBSD_iq'), '-m2'); 
        close

        % calculating grain boundaries and removing defects
        pointsToRemove = (ebsd.iq < (0.12*max(ebsd.iq)));
        ebsd(pointsToRemove) = 'notIndexed';

        [grains, ebsd.grainId] = calcGrains(ebsd,'angle', 10*degree, 'boundary', 'tight');
        ebsd(grains(grains.grainSize < 3)) = [];
        F = halfQuadraticFilter; % denoise
        F.alpha = 0.25;
        ebsd = smooth(ebsd,F,'fill',grains);
        [grains, ebsd.grainId, ebsd.mis2mean] = calcGrains(ebsd,'angle', 10*degree, 'boundary', 'tight');
        grains = smooth(grains, 5);

        % plotting EBSD w/grain boundaries
        figure;
        plot(ebsd,ebsd.orientations)
        hold on
        plot(grains.boundary, 'lineWidth', 1)
        export_fig(sprintf('%s %s.png', samp_name, ' EBSD wGBs'), '-m2'); 
        close

        % saving processed variables for later
        save(samp_name,'ebsd','grains');
    end
end

% NEEDED >> Add way to select specific sets of data and load in their .mat files

%% Tony initial version - create an EBSD variable containing the data
% ebsd = EBSD.load(fname,csFCC,'interface','ang', 'convertEuler2SpatialReferenceFrame','setting 2'); 
% ebsd = rotate(ebsd,180*degree,'keepEuler'); % rotate to have build direction aligned y-axis
% ebsd_raw = ebsd;
% figure
% plot(ebsd_raw,ebsd_raw.orientations,'micronbar','off')
% figure
% plot(ebsd_raw,ebsd_raw.prop.iq,'micronbar','off')
% colormap('gray')
% 
% [grains,ebsd.grainId] = calcGrains(ebsd,'angle',10*degree, 'minPixel',50);
% F = halfQuadraticFilter; % denoise
% F.alpha = 0.25;
% ebsd = smooth(ebsd,F,'fill',grains);
% % remove small grains and make them unindexed
% toRemove = grains.grainSize<15; 
% ebsd(grains(toRemove)) = 'notIndexed';
% [grains, ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'angle', 10*degree, 'minPixel',50);
% 
% % some of false grains appear inside unindexed region
% % better to avoid for loop later
% falsegrains = [];
% for z = 1:length(grains)
%     tempId = grains(z).neighbors('full');
%     temp = height(tempId) == 1 && any(grains(tempId).isIndexed==0);
%     falsegrains(z) = temp;
% end
% ebsd(grains(falsegrains'==1)) = 'notIndexed';
% [grains, ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'angle', 10*degree);
% figure
% plot(ebsd,ebsd.orientations,'micronbar','off')
% figure
% plot(grains,'micronbar','off')
% 
% figure
% plot(ebsd,ebsd.orientations,'micronbar','off')
% figure
% plot(grains,'micronbar','off')
% 
% % save processed data. once have, no need to run code above
% save(outname,'ebsd','grains');

%% Data analysis
% [grains_for_PP, ebsd.grainId,ebsd2.mis2mean] = calcGrains(ebsd,'unitCell'); % strategically have unindexed as phase for post data processing
% [grains, ebsd.grainId,ebsd2.mis2mean] = calcGrains(ebsd('indexed'),'unitCell');

% Grain orientation spread. 
gos = grains.GOS./degree;
figure
plot(grains,gos,'micronbar','off')
setColorRange([0 8])

fraction = 0.4;
numPoints = round(fraction*length(gos));
Ran_GOS_id = randperm(length(gos),numPoints);
dilutedGOS = gos(Ran_GOS_id);

figure
h = daviolinplot(dilutedGOS,'box',0,...
'boxcolor','w','scatter',2,'jitter',1,'scattercolor',[0.98, 0.40, 0.35],'colors',[0.98, 0.40, 0.35],...
'scattersize',2,'scatteralpha',0.7,'linkline',1,'outliers',0,...
'xtlabels', ' '); 
ylabel('GOS (\circ)');
xlabel(' ')
ylim([0 8])
xl = xlim; xlim([xl(1)-0.1, xl(2)+0.4]); % make more space for the legend
set(h.sc,'MarkerEdgeColor','none');      % remove marker edge color
set(gca,'FontSize',12);
set(gca,'Box','on')
set(gca, 'LineWidth', 2)

% merge twins for fcc
twinning = orientation('axis', Miller(1,1,1,ebsd.CS),'angle',60*degree,ebsd.CS)\orientation('euler',0,0,0,ebsd.CS);
gB = grains.boundary;
isTwinning = angle(gB.misorientation,twinning) < 8.86*degree;
twinBoundary = gB(isTwinning);
[mergedGrains,parentId] = merge(grains('indexed'),ebsd,twinning);
figure
plot(grains,grains.meanOrientation)
hold on
plot(twinBoundary,'lineWidth',2,'lineColor','w')
hold off
% exclude grains not fully inside scanned region
% outerBoundaries_id = any(mergedGrains.boundary.grainId==0,2); % exclude grains not fully inside scanned region
% grains_exclude_id = mergedGrains.boundary(outerBoundaries_id).grainId;
% grains_exclude_id(grains_exclude_id==0) = [];
% mergedGrains(grains_exclude_id) = [];

% Equivalent grain diameter
grain_diameter = 2.*equivalentRadius(mergedGrains);
figure
plot(mergedGrains,mergedGrains.meanOrientation,'micronbar','off')
figure
boxplot(grain_diameter)
ylabel('Grain equivalent diameter \mum')
PrettyPlotsSingle

% Grain shape
ShapeFactor = grains.shapeFactor;
figure
plot(grains,ShapeFactor)
setColorRange([0 3])

%% Pole figures for texture analysis
CS = ebsd.CS;
odf = calcDensity(ebsd('Face Centered Cubic').orientations);

% change direction to have build direction upward later
% plotx2east;
% plotzIntoPlane;

setMTEXpref('FontSize',40)
pfAnnotations = @(varargin) text([vector3d.X,-vector3d.Y], {'', ''},...% put RD ND TD you want for X and Y
    'BackgroundColor', 'w', 'FontSize',32, 'tag', 'axesLabels', varargin{:});
setMTEXpref('pfAnnotations', pfAnnotations);


pt = 2;
pos_name = {'ND_Offset', 'Center', 'TD_Offset', 'RD_Offset'};

figure
plotPDF(odf, Miller({1,0,0},CS),'antipodal','projection','eangle');
plotPDF(odf, Miller({1,0,0},CS),'antipodal', 'contour', 0.3:0.3:4, 'linewidth', 1.5, 'linecolor', 'black', 'ShowText', false, 'labelspacing', 300, 'projection','eangle', 'add2all');
mtexColorMap parula
setColorRange([0 2])
mtexColorbar

figure
plotPDF(odf, Miller({1,1,0},CS),'antipodal','projection','eangle');
plotPDF(odf, Miller({1,1,0},CS),'antipodal', 'contour', 0.3:0.3:4, 'linewidth', 1.5, 'linecolor', 'black', 'ShowText', false, 'labelspacing', 300, 'projection','eangle', 'add2all');
mtexColorMap parula
setColorRange([0 2])
mtexColorbar

figure
plotPDF(odf, Miller({1,1,1},CS),'antipodal','projection','eangle');
plotPDF(odf, Miller({1,1,1},CS),'antipodal', 'contour', 0.3:0.3:4, 'linewidth', 1.5, 'linecolor', 'black', 'ShowText', false, 'labelspacing', 300, 'projection','eangle', 'add2all');
mtexColorMap parula
setColorRange([0 2])
mtexColorbar

%% Grains near defects
defect_id = grains_for_PP.isIndexed==0 & grains_for_PP.grainSize > 100;
defect_grains = grains_for_PP(defect_id);
figure
plot(grains_for_PP,grains_for_PP.meanOrientation,'micronbar','off')
hold on
plot(defect_grains.boundary,'linewidth',2)
<<<<<<< HEAD
hold off
=======
hold off
>>>>>>> 9b347aabc59a987e1c3fdb59b807f1d6b922f579
