sourceFolder = '/home/erfan/AD_UCSF/Results/source/';

% Create a new figure to hold all the subplots
figure('Position', [10, 10, 1200, 1200]);

% Get a list of all directories in the source folder
dirList = dir(sourceFolder);
dirList = dirList([dirList.isdir]);
dirList = dirList(~ismember({dirList.name}, {'.', '..'}));

% Preallocate the matrix to store all the source ROI data
filt.low = [5];
filt.high = [50];
fs = 200;

% Calculate the number of subplots needed and the grid size
numPlots = numel(dirList);
numRows = ceil(sqrt(numPlots)); % Number of rows
numCols = ceil(numPlots / numRows); % Number of columns

% Loop through each directory
for i = 1:numPlots
    % Get the current directory path
    currentDir = fullfile(sourceFolder, dirList(i).name);

    % Load the source_roi_data field
    load(fullfile(currentDir, 'source_rec_results.mat'), 'source_roi_data');

    % Process the data
    source_roi_data = source_roi_data(:, :, :);
    [b_orig, b_anti, b_orig_norm, b_anti_norm] = fp_pac_bispec(source_roi_data, fs, filt);

    % Create a subplot in the calculated grid format
    subplot(numRows, numCols, i);
    imagesc(b_anti);
    % xticks=labels;
%     clim([0 1]);
    colormap hot;
    title(dirList(i).name);
%     set(gca, 'XTick', 1:numel(labels), 'XTickLabel', labels, 'XTickLabelRotation', 45);  % Rotate if needed
%     set(gca, 'YTick', 1:numel(labels), 'YTickLabel', labels);

end

% Save the figure with all subplots
saveas(gcf, '/home/erfan/AD_UCSF/Results/all_figures.png');
