
for idx = 1:length(matNames)
    selectedMatFile = matNames{idx};
    load(selectedMatFile)
    
    for i = 1:length(allFileInfo)
        if strcmpi(allFileInfo{i}.sampleName, allMatBases{idx})
            sampleName = allFileInfo{i}.sampleName;
            thisFile = allFileInfo{i};
            save(strcat(sampleName, '_wSubGrains'),'ebsd','grains','thisFile');
        end
    end
    
    copy_ebsd = ebsd;
    copy_grains = grains;
    
    subGB = grains.innerBoundary;
    connectedGB = subGB(subGB.componentSize > 3);
    
    % leftIQ  = ebsd(connectedGB.ebsdId(:,1)).iq;
    % rightIQ = ebsd(connectedGB.ebsdId(:,2)).iq;
    % avgIQ   = (leftIQ + rightIQ) / 2;
    % connectedGB = connectedGB(avgIQ > 0.4 * max(ebsd.iq));
    
    figure;
    plot(ebsd,ebsd.orientations)
    hold on
    plot(grains.boundary, 'lineWidth', 1)
    hold on
    plot(connectedGB, 'lineWidth', 0.75, 'lineColor','blue')
    export_fig(char(fullfile(imgPath, sprintf('%s %s.png', sampleName, ' EBSD wSubGBs'))), '-m2'); 
    close
    
    % GROD
    grod = ebsd.calcGROD(grains);
    
    % plot the misorientation angle of the GROD
    figure;
    plot(ebsd,grod.angle./degree,'micronbar','off')
    setColorRange([0,20])
    mtexColorbar('title',{'misorientation angle in degree'})
    mtexColorMap LaboTeX
    
    % overlay grain and sub-grain boundaries
    hold on
    plot(grains.boundary,'lineWidth',1)
    plot(grains.innerBoundary,'lineWidth',0.75,'edgeAlpha',grains.innerBoundary.misorientation.angle / (5*degree))
    hold off
    export_fig(char(fullfile(imgPath, sprintf('%s %s.png', sampleName, ' GROD'))), '-m2'); 
    close
    
    % KAM
    ebsd_grid = ebsd.gridify;
    kam = ebsd_grid.KAM / degree;
    
    figure;
    plot(ebsd_grid,kam,'micronbar','off')
    setColorRange([0,5])
    mtexColorbar('title',{'degree'})
    mtexColorMap LaboTeX
    hold on
    plot(grains.boundary,'lineWidth',1)
    hold off
    export_fig(char(fullfile(imgPath, sprintf('%s %s.png', sampleName, ' KAM_5deg'))), '-m2');
    close
    
    figure;
    plot(ebsd_grid,ebsd_grid.KAM('threshold',1*degree) ./ degree,'micronbar','off')
    setColorRange([0,1])
    mtexColorbar('title',{'degree'})
    mtexColorMap LaboTeX
    hold on
    plot(grains.boundary,'lineWidth',1)
    hold off
    export_fig(char(fullfile(imgPath, sprintf('%s %s.png', sampleName, ' KAM_thres1deg'))), '-m2');
    close
    
    % Planning out code for different grain/subgrain analysis metrics
    gb = grains.boundary('fcc', 'fcc');
    
    xMax = max(ebsd.prop.x); % width of ebsd scan (x)
    yMax = max(ebsd.prop.y); % height of ebsd scan (y)
    ebsd_Area = xMax*yMax; % area um^2
    
    grain_ED = 2*grains.equivalentRadius;
    subGBDensity_ebsdArea = sum(connectedGB.segLength)/ebsd_Area; % units: um^-1
    subGBLengthPerGrain = accumarray(connectedGB.grainId(:), repmat(connectedGB.segLength, 2, 1), [length(grains) 1]);
    subGBLengthPerGrainSize = accumarray(connectedGB.grainId(:), repmat(connectedGB.segLength, 2, 1), [length(grain_ED) 1]);
    
    % removing edge grains from consideration
    outerBoundaries_id = any(grains.boundary.grainId==0,2);
    grain_id = grains.boundary(outerBoundaries_id).grainId;
    grain_id(grain_id==0) = [];
    grains(grain_id) = [];
    
    isFCCGrain = grains.phaseId == 2;
    subGBLengthPerGrain = subGBLengthPerGrain(isFCCGrain);
    subGBLengthPerGrainSize = subGBLengthPerGrainSize(isFCCGrain);
    
    grains_FCC = grains(isFCCGrain);
    grains_notInd = grains(~isFCCGrain);
    
    grains_FCC_ED = grain_ED(isFCCGrain);
    grains_FCC_AR = grains_FCC.aspectRatio;
    
    grains_notInd(grains_notInd.grainSize < 150);
    grains_notInd_ED = 2*grains_notInd.equivalentRadius;
    grains_notInd_AR = grains_notInd.aspectRatio;
    
    % GOS
    gos = grains_FCC.GOS./degree;
    figure;
    plot(grains_FCC,gos,'micronbar','off')
    setColorRange([0 10])
    mtexColorbar('title',{'degree'})
    export_fig(char(fullfile(imgPath, sprintf('%s %s.png', sampleName, ' GOS'))), '-m2'); 
    close
    
    % sphericity
    Psi = grains_FCC.area ./ grains_FCC.perimeter('withInclusion') ./ grains_FCC.equivalentRadius;
    
    figure;
    plot(grains_FCC, Psi, 'colorrange', [0 0.5],'micronbar','off')
    
    mtexColorbar ('title','Sphericity Parameter')
    mtexColorMap jet
    export_fig(char(fullfile(imgPath, sprintf('%s %s.png', sampleName, ' sphericity'))), '-m2'); 
    close
    
    % Saving results
    GrainResults.subGrains = connectedGB;
    GrainResults.subGBDensity_ebsdArea = subGBDensity_ebsdArea;
    GrainResults.subGBLengthPerGrain = subGBLengthPerGrain;
    GrainResults.subGBLengthPerGrainSize = subGBLengthPerGrainSize;
    GrainResults.equDia_FCC = grains_FCC_ED;
    GrainResults.aspRatio_FCC = grains_FCC_AR;
    GrainResults.sphericity_FCC = Psi;
    GrainResults.equDia_notIndx = grains_notInd_ED;
    GrainResults.aspRatio_notIndx = grains_notInd_AR;
    GrainResults.GROD_angle = grod.angle./degree;
    GrainResults.KAM = ebsd_grid.KAM;
    GrainResults.GOS = gos;
    
    ebsd = copy_ebsd;
    grains = copy_grains;
    save(strcat(sampleName, '_wSubGrains'),'ebsd','grains','thisFile', 'GrainResults');
end