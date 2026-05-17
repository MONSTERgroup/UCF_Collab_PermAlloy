%% UCF Collab EBSD Analysis

% Add path to Monster Utilities
addpath 'C:\Users\leey2\Documents\GitHub\EBSD-EDS-combination'
addpath 'C:\Users\leey2\Documents\HeRX-utilities'

% crystal symmetry
CS = {... 
  'notIndexed',...
  crystalSymmetry('432', [4 4 4], 'mineral', 'Face Centered Cubic', 'color', [0.53 0.81 0.98])};

prefix = 'C:\Users\leey2\University of Florida\MONSTER group - Documents\Lab Data\AM PermAlloy (UCF Collab)\Test_run';  

pname_subroot = {'400',...
    '900',...
    '1200'};

pname_sub_subroot = {'XY',...
    'XZ'};

% loop over deformation levels
for ii = 3:3
    for tt = 1:1

        pname = [prefix filesep '05.13.26 PermAlloy As-Built ' pname_subroot{ii} 'mms 200W ' pname_sub_subroot{tt}]; 
        
        % loop over quads
        for jj = 1:1
        
        
            %% Specify File Names
            
            clear ebsd ebsd2 grains data paddedData filteredData Grain_Avg_EDS Rollback toRemove toFlip
            
            % which files to be imported
            fname = [pname ' Scan ' num2str(jj) '.ang']; % later save name differently, _ instead of space, and make more concisely
            
            % outname = [pname filesep 'Scan_' num2str(jj) 'new' '.mat'];
            %% Import the Data
            
            % create an EBSD variable containing the data
            ebsd = EBSD.load(fname,CS,'interface','ang',...
              'convertEuler2SpatialReferenceFrame','setting 3'); % setting 3 to rotate 90
            ebsd_raw = ebsd;
            figure
            plot(ebsd_raw,ebsd_raw.orientations,'micronbar','off')

            [grains,ebsd.grainId] = calcGrains(ebsd,'angle',10*degree, 'minPixel',15);
            F = halfQuadraticFilter; % denoise
            F.alpha = 0.25;
            ebsd = smooth(ebsd,F,'fill',grains);
            % remove small grains and make them unindexed
            toRemove = grains.grainSize<15; 
            ebsd(grains(toRemove)) = 'notIndexed';
            [grains, ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'angle', 10*degree, 'minPixel',15);

            % some of false grains appear inside unindexed region
            % better to avoid for loop later
            falsegrains = [];
            for z = 1:length(grains)
                tempId = grains(z).neighbors('full');
                temp = height(tempId) == 1 && any(grains(tempId).isIndexed==0);
                falsegrains(z) = temp;
            end
            ebsd(grains(falsegrains'==1)) = 'notIndexed';
            [grains, ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'angle', 10*degree);
            figure
            plot(ebsd_raw,ebsd_raw.orientations,'micronbar','off')
            figure
            plot(grains,'micronbar','off')

            % small points on scratches, filling them
            toFlip = grains.grainSize<10 & grains.isIndexed==0; 
            ebsd(grains(toFlip)).phaseId = 2;
            ebsd(grains(toFlip)).phase = 0;
            [grains, ebsd.grainId,ebsd2.mis2mean] = calcGrains(ebsd,'angle', 10*degree, 'minPixel',15);
            ebsd = smooth(ebsd,grains);
            ebsd('notIndexed') = [];
            
            figure
            plot(ebsd_raw,ebsd_raw.orientations,'micronbar','off')
            figure
            plot(grains,'micronbar','off')
    
        end
    end
end