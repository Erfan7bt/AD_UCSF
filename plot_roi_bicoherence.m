% Define the directory containing the results
resultsDir = '/home/erfan/AD_UCSF/Results/source/';
save_dir = '/home/erfan/AD_UCSF/roi_pac/';

% Get a list of all subjects in the directory
subjectFiles = dir(fullfile(resultsDir));
% Remove the '.' and '..' entries and files
subjectFiles = subjectFiles(~ismember({subjectFiles.name}, {'.', '..'}));

% Loop through each subject and load the source_rec_results
for i = 1:length(subjectFiles)
    % Construct the full path to the subject's results file
    subjectFilePath = fullfile(resultsDir, subjectFiles(i).name, 'source_rec_results.mat');
    
    load(subjectFilePath, 'source_roi_data', 'labels');

    T = 2;
    fs = 600;
    epleng = T * fs; % 2 sec
    for segleng = [epleng epleng/2]
        for segshift = [segleng segleng/2]
            df = 1 / (segleng / fs);
            desired_f = 50;

            maxfreqbins = desired_f / df + 1;    
            for j = 1%:size(source_roi_data, 1)
                b= {};
                data = source_roi_data(j, :, :);
                data = data(:,:)';
                [cs, csnr, nave] = data2bs_univar(data, segleng, segshift, epleng, maxfreqbins);
                
                cs = abs(cs);
                b{1}= reshape(cs(1, :, :), size(cs, 2), []);
                cs = cs ./ csnr;
                b{2}= reshape(cs(1, :, :), size(cs, 2), []);

                for B = 1:2
                    fig = figure('Visible', 'on');
                    imagesc(b{B});
                    colormap hot
                    colorbar;
                    name = labels{j};
                    name = [name ' segl' num2str(segleng) ' segs' num2str(segshift)];
                    title(name);
                    xticks = linspace(1, maxfreqbins, 6);  % 6 evenly spaced ticks (0, 10, 20, 30, 40, 50 Hz)
                    xticklabels = 0:10:desired_f;  % Corresponding labels from 0 to 50 Hz
                    % Apply the xticks and corresponding labels
                    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
                    set(gca, 'YTick', xticks, 'YTickLabel', xticklabels);  % Same for y-axis if needed
                    % Add labels and title (optional)
                    xlabel('Frequency (Hz)');
                    ylabel('Frequency (Hz)');

                    if ~isfolder([save_dir subjectFiles(i).name])
                        mkdir([save_dir subjectFiles(i).name]);
                        mkdir([save_dir subjectFiles(i).name '/Bicoh']);
                    end
                    if B == 2
                        clim([0 1])
                        saveas(fig, [save_dir subjectFiles(i).name '/Bicoh/' name '.png']);
                        disp(['Saved ' name '.png']);
                        close(fig);
                    else
                        saveas(fig, [save_dir subjectFiles(i).name '/' name '.png']);
                        disp(['Saved ' name '.png']);
                        close(fig);
                    end
                end
            end
        end
    end
 end
