% Demo: Sensitivity Analysis of Grain Metrics to Angular Threshold
% Using Gaussian Process Regression

%% 1. Generate Synthetic Dataset (replace with your real data loading)
thresholds = 5:2:25;  % Angular thresholds: 5, 7, 9, ..., 25 degrees
n_samples = 8;        % Number of EBSD maps/samples
n_thresholds = length(thresholds);

% Pre-allocate results
data = [];

for sample_id = 1:n_samples
    for i = 1:n_thresholds
        thresh = thresholds(i);
        
        % Simulate grain reconstruction metrics (replace with your real analysis)
        metrics = simulateGrainMetrics(thresh, sample_id);
        
        % Store results
        row = [sample_id, thresh, metrics.mean_grain_size, ...
               metrics.gos_mean, metrics.gos_std, metrics.n_grains];
        data = [data; row];
    end
end

% Convert to table for easier handling
colNames = {'SampleID', 'Threshold', 'MeanGrainSize', 'GOS_Mean', 'GOS_Std', 'NumGrains'};
T = array2table(data, 'VariableNames', colNames);

%% 2. Fit Gaussian Process for each metric
metrics_to_analyze = {'MeanGrainSize', 'GOS_Mean', 'GOS_Std', 'NumGrains'};
gp_models = containers.Map();
results = struct();

% Fine grid for prediction
thresh_fine = linspace(min(thresholds), max(thresholds), 100);

for m = 1:length(metrics_to_analyze)
    metric_name = metrics_to_analyze{m};
    
    % Prepare data (average across samples for simplicity)
    metric_data = grpstats(T, 'Threshold', 'mean', 'DataVars', metric_name);
    X = metric_data.Threshold;
    y = metric_data.(['mean_' metric_name]);
    
    % Fit GP (requires Statistics and Machine Learning Toolbox)
    gpr = fitrgp(X, y, 'KernelFunction', 'squaredexponential', ...
                 'OptimizeHyperparameters', 'auto', 'Verbose', 0);
    
    % Predictions with uncertainty
    [y_pred, y_std] = predict(gpr, thresh_fine');
    
    % Store results
    results.(metric_name).model = gpr;
    results.(metric_name).thresh_fine = thresh_fine;
    results.(metric_name).y_pred = y_pred;
    results.(metric_name).y_std = y_std;
    results.(metric_name).raw_data = [X, y];
end

%% 3. Compute Sensitivity Metrics
sensitivity_report = struct();

for m = 1:length(metrics_to_analyze)
    metric_name = metrics_to_analyze{m};
    
    % Get predictions
    thresh = results.(metric_name).thresh_fine;
    y_pred = results.(metric_name).y_pred;
    y_std = results.(metric_name).y_std;
    
    % Compute derivative (sensitivity per degree)
    dy_dt = gradient(y_pred, thresh);
    
    % Key sensitivity metrics
    baseline_val = y_pred(1);  % Value at minimum threshold
    change_10_to_15 = interp1(thresh, y_pred, 15) - interp1(thresh, y_pred, 10);
    relative_change = change_10_to_15 / interp1(thresh, y_pred, 10);
    max_sensitivity = max(abs(dy_dt));
    
    % Find "stable region" (where sensitivity is low)
    stability_threshold = 0.1 * max_sensitivity;
    stable_indices = abs(dy_dt) < stability_threshold;
    
    % Store sensitivity metrics
    sensitivity_report.(metric_name).change_10_to_15 = change_10_to_15;
    sensitivity_report.(metric_name).relative_change = relative_change;
    sensitivity_report.(metric_name).max_sensitivity = max_sensitivity;
    sensitivity_report.(metric_name).stable_range = [min(thresh(stable_indices)), ...
                                                     max(thresh(stable_indices))];
end

%% 4. Visualization
figure('Position', [100, 100, 1200, 800]);

for m = 1:length(metrics_to_analyze)
    subplot(2, 2, m);
    metric_name = metrics_to_analyze{m};
    
    % Get data
    thresh = results.(metric_name).thresh_fine;
    y_pred = results.(metric_name).y_pred;
    y_std = results.(metric_name).y_std;
    raw_data = results.(metric_name).raw_data;
    
    % Plot GP prediction with uncertainty
    fill([thresh, fliplr(thresh)], [y_pred + 1.96*y_std; flipud(y_pred - 1.96*y_std)], ...
         [0.7 0.7 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    hold on;
    plot(thresh, y_pred, 'b-', 'LineWidth', 2);
    scatter(raw_data(:,1), raw_data(:,2), 50, 'ro', 'filled');
    
    % Highlight 10° and 15° points
    plot([10 10], ylim, 'g--', 'LineWidth', 1.5);
    plot([15 15], ylim, 'r--', 'LineWidth', 1.5);
    
    xlabel('Angular Threshold (degrees)');
    ylabel(strrep(metric_name, '_', ' '));
    title(sprintf('%s vs Threshold', strrep(metric_name, '_', ' ')));
    grid on;
    legend('95% CI', 'GP Prediction', 'Data Points', '10°', '15°', 'Location', 'best');
end

sgtitle('Grain Metrics Sensitivity to Angular Threshold');

%% 5. Print Sensitivity Report
fprintf('\n=== SENSITIVITY ANALYSIS REPORT ===\n\n');
for m = 1:length(metrics_to_analyze)
    metric_name = metrics_to_analyze{m};
    rep = sensitivity_report.(metric_name);
    
    fprintf('Metric: %s\n', strrep(metric_name, '_', ' '));
    fprintf('  Change (10° → 15°): %.3f\n', rep.change_10_to_15);
    fprintf('  Relative change: %.1f%%\n', rep.relative_change * 100);
    fprintf('  Max sensitivity: %.3f per degree\n', rep.max_sensitivity);
    fprintf('  Stable range: %.1f° - %.1f°\n\n', rep.stable_range);
end

%% Helper Function: Simulate grain metrics (replace with your real analysis)
function metrics = simulateGrainMetrics(threshold, sample_id)
    % Simulate realistic grain reconstruction behavior
    % Replace this entire function with your actual grain analysis
    
    % Add some sample-to-sample variation
    rng(sample_id + 100);  % Reproducible but different per sample
    
    % Realistic trends: higher threshold → larger grains, fewer grains
    base_grain_size = 15 + 3 * (threshold/10) + 2*randn();
    base_n_grains = max(50, 200 - 5*threshold + 10*randn());
    
    % GOS tends to decrease with higher threshold (less fragmentation)
    base_gos_mean = max(0.5, 3 - 0.1*threshold + 0.3*randn());
    base_gos_std = max(0.2, 1.5 - 0.05*threshold + 0.2*randn());
    
    metrics.mean_grain_size = max(1, base_grain_size);
    metrics.n_grains = round(max(10, base_n_grains));
    metrics.gos_mean = max(0.1, base_gos_mean);
    metrics.gos_std = max(0.05, base_gos_std);
end