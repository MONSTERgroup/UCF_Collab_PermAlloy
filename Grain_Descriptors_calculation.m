%% Grain Boundary Threshold. Check Descriptors' sensitivity
% Get list of all .mat files in current directory

folder = "C:\Users\leey2\University of Florida\MONSTER group - Documents\Lab Data\AM PermAlloy (UCF Collab)\Test_run\x150 mat files";
mat_files  = dir(fullfile(folder, "*.mat"));

% Initialize variable to store all analyses
angle_thresh = 5:2:15;
Variance_of_dist = {};
PC1 = [];
Grain_GOS_var = [];
Grain_GOS_Mean = [];
Grain_GOS_Median = [];
Grain_Size_var = [];
Grain_Size_Mean = [];
Grain_Size_Median = [];
Grains_ShapeFactor_var = [];
Grains_ShapeFactor_Mean = [];
Grains_ShapeFactor_Median = [];
Grains_Subgrain_var = [];
Grains_Subgrain_Mean = [];
Grains_Subgrain_Median = [];
Grain_Number = [];
% Loop through each .mat file
for k = 6:6%numel(mat_files)
    fpath = fullfile(mat_files(k).folder, mat_files(k).name);

    load(fpath); % .mat file contains grains and ebsd. grains is when thresh angle was 10 deg

    for i = 1:numel(angle_thresh)
        [grains, ebsd.grainId, ebsd.mis2mean] = calcGrains(ebsd,'angle', angle_thresh(i)*degree, 'boundary', 'tight');
        grains = smooth(grains, 5);
        Grain_GOS = grains('Face Centered Cubic').GOS./degree;
        Grain_GOS_Mean = [Grain_GOS_Mean, mean(Grain_GOS)];
        Grain_GOS_Median = [Grain_GOS_Median, median(Grain_GOS)];
        Grain_GOS_var = [Grain_GOS_var, var(Grain_GOS)];
        Grain_Size = 2*equivalentRadius(grains('Face Centered Cubic'));
        Grain_Size_Mean = [Grain_Size_Mean, mean(Grain_Size)];
        Grain_Size_Median = [Grain_Size_Median, median(Grain_Size)];
        Grain_Size_var = [Grain_Size_var, var(Grain_Size)];
        Grains_ShapeFactor = grains('Face Centered Cubic').shapeFactor;
        Grains_ShapeFactor_Mean = [Grains_ShapeFactor_Mean, mean(Grains_ShapeFactor)];
        Grains_ShapeFactor_Median = [Grains_ShapeFactor_Median, median(Grains_ShapeFactor)];
        Grains_ShapeFactor_var = [Grains_ShapeFactor_var, var(Grains_ShapeFactor)];
        Grains_Subgrain = grains('Face Centered Cubic').subBoundaryLength; 
        Grains_Subgrain_Mean = [Grains_Subgrain_Mean, mean(Grains_Subgrain)];
        Grains_Subgrain_Median = [Grains_Subgrain_Median, median(Grains_Subgrain)];
        Grains_Subgrain_var = [Grains_Subgrain_var, var(Grains_Subgrain)];
        Grain_Number = [Grain_Number, length(grains('Face Centered Cubic'))];
        PCA_Matrix = [Grain_GOS, Grain_Size, Grains_ShapeFactor, Grains_Subgrain];
        [coeff,score,latent] = pca(PCA_Matrix); 
        PC1 = [PC1,coeff(:,1)];
        Variance_of_dist{i} = latent;
        clear grains PCA_Matrix
    end
end

%% Measured points at discrete GB thresh, predict at other values GB thresh  
Grain_Number_Gaussian = Grain_Number';               

% Fit model 
gpr = fitrgp(angle_thresh', Grain_Number_Gaussian, ...
             'KernelFunction','squaredexponential', ...
             'Standardize',true);

% Predict on finer grid
angle_thresh_fine = (5:0.1:15)';           % Xnew column
[y_pred, y_sd] = predict(gpr, angle_thresh_fine);

figure;
fill([angle_thresh_fine; flipud(angle_thresh_fine)], ...
[y_pred + 1.96*y_sd; flipud(y_pred - 1.96*y_sd)], ... % 1.96 for CI 95%
[0.7 0.7 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
hold on;
plot(angle_thresh_fine, y_pred, 'b-', 'LineWidth', 2);
scatter(angle_thresh, Grain_Number_Gaussian, 50, 'ro', 'filled');

% GB angle thresh where Descriptors least likely be changing
dydx = gradient(y_pred) ./ gradient(angle_thresh_fine);
d2ydx2 = gradient(dydx) ./ gradient(angle_thresh_fine);
figure
plot(angle_thresh_fine,abs(dydx))
% Find where the 2nd derivative crosses zero
crossings = find(diff(sign(d2ydx2)));

% Interpolate for the exact X & Y coordinates of the inflection point
infl_x = interp1(d2ydx2(crossings), angle_thresh_fine(crossings), 0);

% GB thresh (might be better to be an integer)
GB_Thresh = round(infl_x);