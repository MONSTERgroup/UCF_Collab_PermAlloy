function nDashed = applyDashedLinesToBoxplot(chosenGroup, numFiles, ax)
% Apply dashed line styles to XZ build direction
%
%   nDashed = applyDashedLinesToBoxplot(chosenGroup, numFiles)
%   nDashed = applyDashedLinesToBoxplot(chosenGroup, numFiles, ax)
%
% Inputs:
%   chosenGroup - String: 'allXY', 'allXZ', 'onlyAsBuilt', 'onlyHIP', etc.
%   numFiles    - Integer: Total number of files/boxes in the plot
%   ax          - (Optional) Axes handle. Default: gca
%
% Output:
%   nDashed     - Number of boxes that were made dashed (0 if none)
%
% Description:
%   For groups with mixed build directions (XY and XZ), applies dashed
%   line styles to the second half (XZ direction) to visually distinguish
%   them from the first half (XY direction).
%
% Example:
%   lognormalBoxplot(data, groups);
%   PrettyPlotsSingle_BoxPlot(colors);
%   applyDashedLinesToBoxplot('onlyAsBuilt', 6);

% Handle optional axes argument
if nargin < 3
    ax = gca;
end

% Determine which files need dashed lines based on group
applyDashing = true;
switch chosenGroup
    case 'allXY'
        applyDashing = false;
        nDashed = 0;

    case 'allXZ'
        nDashed = 5;

    case 'onlyAsBuilt'
        nDashed = 3;
        
    case 'onlyHIP'
        nDashed = 2;
        
    case 'compTopVsCenter'
        % Custom handling for top vs center comparison
        nDashed = floor(numFiles / 2);
        
    otherwise
        % Generic: second half dashed
        nDashed = floor(numFiles / 2);
end

if ~applyDashing
    return;
end

allBoxes = findobj(ax, 'Tag', 'Box');
allMedians = findobj(ax, 'Tag', 'Median');
allUpperWhiskers = findobj(ax, 'Tag', 'Upper Whisker');
allLowerWhiskers = findobj(ax, 'Tag', 'Lower Whisker');
allUpperAdjValues = findobj(ax, 'Tag', 'Upper Adjacent Value');
allLowerAdjValues = findobj(ax, 'Tag', 'Lower Adjacent Value');

if nDashed == 0
    return;
end

% Select elements to make dashed (first nDashed in reversed array)
dashedBoxes = allBoxes(1:nDashed);
dashedMedians = allMedians(1:min(nDashed, length(allMedians)));
dashedUpperWhiskers = allUpperWhiskers(1:min(nDashed, length(allUpperWhiskers)));
dashedLowerWhiskers = allLowerWhiskers(1:min(nDashed, length(allLowerWhiskers)));
dashedUpperAdjValues = allUpperAdjValues(1:min(nDashed, length(allUpperAdjValues)));
dashedLowerAdjValues = allLowerAdjValues(1:min(nDashed, length(allLowerAdjValues)));

% Apply dashed line style
set(dashedBoxes, 'LineStyle', ':', 'LineWidth', 2)
set(dashedMedians, 'LineStyle', ':', 'LineWidth', 2)
set(dashedUpperWhiskers, 'LineStyle', ':')
set(dashedLowerWhiskers, 'LineStyle', ':')
set(dashedUpperAdjValues, 'LineStyle', ':')
set(dashedLowerAdjValues, 'LineStyle', ':')

end