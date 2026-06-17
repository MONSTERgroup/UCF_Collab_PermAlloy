fpath = "C:\Users\leey2\University of Florida\MONSTER group - Documents\Lab Data\AM PermAlloy (UCF Collab)\Test_run\x150 mat files";
fname = {'As-Built 900mms 200W XZ Scan 2.mat', 'HIP 120min 900mms 200W XZ Scan 2.mat', 'HIP 180min 900mms 200W XZ Scan 2.mat'};
plotx2east
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','intoPlane'); % the BD would be Y
setMTEXpref('FontSize',14);
colorScheme = {[250/255 70/255 22/255], [0/255 30/255 165/255], [34/255 146/255 76/255]};

% Initialize cell array to store results
Results_900mms = cell(length(fname), 1);

for i = 1:length(fname)
    fprintf('Processing file %d of %d: %s\n', i, length(fname), fname{i});
    
    % Load the current .mat file
    data = load(fullfile(fpath, fname{i}));
    ebsd = data.ebsd;

    % Rotate to have BD y direction
    rot = rotation.byAxisAngle(zvector,90*degree);
    ebsd = rotate(ebsd,rot,'keepEuler');

    grains = data.grains; 
    CS = ebsd.CS;
    twinning = orientation('axis', Miller(1,1,1,CS),'angle',60*degree,CS)\orientation('euler',0,0,0,CS);

    %% Twin Analysis (just in case we wanna do something with twin for HIPed sample later)
    grainsFCC = grains('Face Centered Cubic');
    gB = grainsFCC.boundary('Face Centered Cubic','Face Centered Cubic');
    isTwinning = angle(gB.misorientation,twinning) < 8.66*degree;
    twinBoundary = gB(isTwinning);
    
    if ~isempty(twinBoundary)
        [mergedGrains,parentId] = merge(grainsFCC,twinBoundary);
        % merge twins for GSD
        outerBoundaries_id = any(mergedGrains.boundary.grainId==0,2);
        grain_id = mergedGrains.boundary(outerBoundaries_id).grainId;
        grain_id(grain_id==0) = [];
        mergedGrains(grain_id) = [];
        gArea = grainsFCC.area;
        twinsPerGrain = zeros(mergedGrains.length,1);
        isTwin = true(grainsFCC.length,1);
        
        for j = 1:mergedGrains.length % Changed loop variable to avoid conflict
            % get child ids
            childId = find(parentId==j);
            % cluster grains of similar orientations
            [fId,~] = calcCluster(grainsFCC(childId).meanOrientation,'maxAngle',...
                15*degree,'method','hierarchical','silent');
            % compute area of each cluster
            clusterArea = accumarray(fId,gArea(childId));
            % label the grains of largest cluster as original grain
            [~,fParent] = max(clusterArea);
            isTwin(childId(fId==fParent)) = false;
            % add number of grains to the array
            twinsPerGrain(j) = length(grainsFCC(childId(isTwin(childId))));
        end
        
        % compute the twinning statistics
        twinAreaFraction = sum(area(grainsFCC(isTwin)))/sum(area(grainsFCC)) * 100;
        numTwins = length(grainsFCC(isTwin));
        twinnedFraction = (sum(twinsPerGrain>0)/length(mergedGrains)) * 100;
        [~,~,twinWidths] = grainsFCC(isTwin).fitEllipse;
    else
        twinAreaFraction = 0;
        twinsPerGrain = zeros(grainsFCC.length,1);
        numTwins = 0;
        twinnedFraction = 0;
        twinWidths = [];
    end
    
    % Store results for this file
    Results_900mms{i}.filename = fname{i};
    Results_900mms{i}.ebsd = ebsd;
    Results_900mms{i}.grains = grains;
    Results_900mms{i}.mergedGrains = mergedGrains;
    Results_900mms{i}.Merged_grains_w_twins = twinsPerGrain(twinsPerGrain > 0);
    Results_900mms{i}.Twin_grain_area_frac = twinAreaFraction;
    Results_900mms{i}.twin_grain_aspect_ratio = grainsFCC(isTwin).aspectRatio;
    Results_900mms{i}.twin_grain_fraction = numel(grainsFCC(isTwin).id)/numel(grainsFCC.id);
    Results_900mms{i}.numTwins = numTwins;
    Results_900mms{i}.twinnedFraction = twinnedFraction;
    Results_900mms{i}.twinWidths = twinWidths;
    
    fprintf('Completed analysis for %s\n', fname{i});
end

% Save all results to a .mat file
save('TwinAnalysisResults.mat', 'Results_900mms');

%% Data analysis after basic calculation above
for i = 1:3
    %% IPDF contour
    figure
    plotIPDF(Results_900mms{i}.ebsd('Face Centered Cubic').orientations,zvector,'contourf') % represent how IPF map is basically. currently it's perpendicular to BD
    setColorRange([0 1.5])

    %% PF
    odf = calcDensity(Results_900mms{i}.ebsd('Face Centered Cubic').orientations);

    setMTEXpref('FontSize',40)
    pfAnnotations = @(varargin) text([vector3d.X,-vector3d.Y], {'TD', 'BD'},...% put RD ND TD you want for X and Y
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
    
    %% GOS
    gos = Results_900mms{i}.grains('Face Centered Cubic').GOS./degree;
    fraction = 0.6;
    numPoints = round(fraction*length(gos)); % dilute number of GOS points
    Ran_GOS_id = randperm(length(gos),numPoints);
    dilutedGOS = gos(Ran_GOS_id);
    
    figure
    h = daviolinplot(dilutedGOS,'box',0,...
    'boxcolor','w','scatter',2,'jitter',1,'scattercolor',colorScheme{i},'colors',colorScheme{i},...
    'scattersize',3,'scatteralpha',0.7,'linkline',1,'outliers',0,...
    'xtlabels', ' '); 
    ylabel('GOS (\circ)');
    xlabel(' ')
    ylim([0 5])
    xl = xlim; xlim([xl(1)-0.1, xl(2)+0.4]); % make more space for the legend
    set(h.sc,'MarkerEdgeColor','none');      % remove marker edge color
    set(gca,'FontSize',12);
    set(gca,'Box','on')
    set(gca, 'LineWidth', 2)
    
    %% Reported RX texture in AM (can be helpful descriptor determining RX % with GOS when GOS alone a bit diffuse)
    % Initial <100> (cube) -> <212> (twin related) or <112> (P-orientation).
    % How they are aligned to the build direction is needed to be understood.
    % Initial texture of AM FCC often shows 100 or 110

    % here in XZ scans, since the build direction is along to the in-plane
    % y-axis, yvector used. The known possible orientation in AM above is
    % || to BD.
    
    threshold = 10; % in degree;
    oriRef = orientation.cube(CS,specimenSymmetry('-1')); 
    
    h100 = [Miller({1,0,0},CS)...
        Miller({-1,0,0},CS)...
        Miller({0,1,0},CS)...
        Miller({0,-1,0},CS)...
        Miller({0,0,1},CS)...
        Miller({0,0,-1},CS)];
    for j = 1:3
        misor = angle(oriRef,Results_900mms{j}.ebsd('Face Centered Cubic').orientations) ./ degree;
        isCube = misor < threshold;
        Results_900mms{j}.fracCube = sum(isCube,'all')/length(isCube(:));
        
        clear r
        for k = 1:length(h100)
            r(k,:) = Results_900mms{j}.ebsd('Face Centered Cubic').orientations * h100(k);
        end
    
        misor2 = min(angle(yvector, r)) ./ degree; % along to the BD so yvector
        is001 = misor2 < threshold; 
        Results_900mms{j}.frac001 = sum(is001)/length(is001); 
    end
    
    u = Miller(1,1,0,CS,'uvw');
    h110 = symmetrise(u);
    h110 = h110';
    for j = 1:3
        clear r
        for k = 1:length(h110)
            r(k,:) = Results_900mms{j}.ebsd('Face Centered Cubic').orientations * h110(k);
        end
    
        misor110 = min(angle(yvector, r)) ./ degree; 
        is011 = misor110 < threshold; 
        Results_900mms{j}.frac011 = sum(is011)/length(is011); 
    end
    
    u = Miller(2,1,2,CS,'uvw'); % one of RX ori
    h212 = symmetrise(u);
    h212 = h212';
    for j = 1:3
        clear r
        for k = 1:length(h212)
            r(k,:) = Results_900mms{j}.ebsd('Face Centered Cubic').orientations * h212(k);
        end
    
        misor212 = min(angle(yvector, r)) ./ degree; 
        is212 = misor212 < threshold; 
        Results_900mms{j}.frac212 = sum(is212)/length(is212); 
    end
    
    u = Miller(1,1,2,CS,'uvw'); % the other RX ori
    h112 = symmetrise(u);
    h112 = h112';
    for j = 1:3
        clear r
        for k = 1:length(h112)
            r(k,:) = Results_900mms{j}.ebsd('Face Centered Cubic').orientations * h112(k);
        end
    
        misor112 = min(angle(yvector, r)) ./ degree; 
        is112 = misor112 < threshold; 
        Results_900mms{j}.frac112 = sum(is112)/length(is112); 
    end
end

%% Grain size distribution
for i = 1:3
    figure
    boxplot(2*Results_900mms{i}.grains.equivalentRadius,'Colors',colorScheme{i})
    ylim([-10 160])
    PrettyPlotsSingle
end