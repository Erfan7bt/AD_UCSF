% Set the main directory path
mainDir = './Figures'; 

% Get the list of all subdirectories
subDirs = dir(mainDir);
subDirs = subDirs([subDirs.isdir]); % Only keep directories
subDirs = subDirs(~ismember({subDirs.name}, {'.', '..'})); % Exclude '.' and '..'

% Initialize a PDF
pdfFileName = 'output_headmodeling.pdf';
append = false; % Set to false initially to create a new PDF

for i = 1:length(subDirs)
    % Get the current subdirectory
    subDirName = subDirs(i).name;
    subDirPath = fullfile(mainDir, subDirName);
    
    % Get all .tif files in the current subdirectory
    tifFiles = dir(fullfile(subDirPath, '*.tif'));
    
    % Group files by their common base name (e.g., "image1_a.tif", "image1_b.tif" -> "image1")
    baseNames = cellfun(@(x) regexp(x, '^[^_]+', 'match', 'once'), {tifFiles.name}, 'UniformOutput', false);
    uniqueBaseNames = unique(baseNames);
    
    for j = 1:length(uniqueBaseNames)
        % Find indices of files that share the same base name
        groupIdx = strcmp(baseNames, uniqueBaseNames{j});
        groupFiles = tifFiles(groupIdx);
        
        % Create a figure for this group of images
        fig = figure('Visible', 'off'); % Create an invisible figure
        
        % Add a text title at the top with the subdirectory name
        annotation('textbox', [0, 0.95, 1, 0.05], 'String', subDirName, ...
                   'FontSize', 10, 'HorizontalAlignment', 'center', 'EdgeColor', 'none');
        
        % Create a tiled layout to place images in a grid
        % Adjust the 'Padding' property to add space between the title and figures
        t = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'loose'); % 'loose' adds space around the grid
        
        % Adjust the position of the layout to add more space
        t.Position(2) = t.Position(2) - 0.05; % Move the tiled layout down by increasing its y-position
        
        % Loop over the grouped files to display them on the same page
        for k = 1:length(groupFiles)
            fileName = groupFiles(k).name;
            filePath = fullfile(subDirPath, fileName);
            
            img = imread(filePath);
            nexttile;
            imshow(img);
            title(fileName, 'Interpreter', 'none'); % Use filename as the caption
        end
        
        % Save the figure content to the PDF
        exportgraphics(fig, pdfFileName, 'Append', append);
        close(fig);
        
        % After the first loop, set append to true
        append = true;
    end
end

disp('PDF created successfully.');
