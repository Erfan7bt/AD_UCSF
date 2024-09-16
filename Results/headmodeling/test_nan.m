% Define the root directory
rootDir = '/home/erfan/AD_UCSF/Results/headmodeling';

% Get a list of all subdirectories
subDirs = dir(rootDir);
subDirs = subDirs([subDirs.isdir]);
subDirs = subDirs(~ismember({subDirs.name}, {'.', '..'}));

% Loop through each subdirectory
for i = 1:numel(subDirs)
    subDirPath = fullfile(rootDir, subDirs(i).name);
    
    % Load the bs_results.mat file
    bsResultsFile = fullfile(subDirPath, 'bs_results.mat');
    load(bsResultsFile, 'leadfield');
    
    % Check the number of NaN values in the leadfield
    numNaN = sum(isnan(leadfield(:)));
    
    % Display the result
    fprintf('Subdirectory: %s\n', subDirPath);
    fprintf('percentage of NaN values: %.2f%%\n', numNaN / numel(leadfield) * 100);
    fprintf('\n');
end