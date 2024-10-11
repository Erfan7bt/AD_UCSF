resultsDir = '/home/erfan/AD_UCSF/Results/source/';
save_dir_main = '/home/erfan/AD_UCSF/hipp_cort/';

% Get a list of all subjects in the directory
subjectFiles = dir(fullfile(resultsDir));
% Remove the '.' and '..' entries and files
subjectFiles = subjectFiles(~ismember({subjectFiles.name}, {'.', '..'}));

% Loop through each subject and load the source_rec_results
for i = 1:length(subjectFiles)
    % Construct the full path to the subject's results file
    subjectFilePath = fullfile(resultsDir, subjectFiles(i).name, 'source_rec_results.mat');

    load(subjectFilePath, 'source_roi_data', 'labels','regions_cortex');
    
    save_dir = [save_dir_main subjectFiles(i).name '/'];

    if ~isfolder(save_dir)
        mkdir(save_dir)
    end
    froi_id = 11; % Hippocampus R
    cortex_id = 16:83;

    fs = 600;          % Sampling frequency
    freq = 1:50-1;     % Frequency range
    numFreqs = length(freq);

    big = zeros(length(cortex_id), numFreqs, numFreqs);

    % Parallel for loop for cortex_id
    for ids_idx = 1:length(cortex_id)
        ids = cortex_id(ids_idx);
        disp(ids);
        tic;

        % Load the data for froi_id and the current cortex id
        data = source_roi_data([froi_id, ids], :, :);

        % Parallelize the frequency computations
        for low_freqs = freq
            for high_freqs = low_freqs+1:numFreqs
                filt.low = low_freqs;
                filt.high = high_freqs;
                %             disp(filt)
                % Call the fp_pac_bispec function
                [b_orig, b_anti, b_orig_norm, b_anti_norm] = fp_pac_bispec(data, fs, filt);

                % Update the 'big' matrix with the computed values
                big(ids_idx, low_freqs, high_freqs) = b_anti_norm(1, 2);
                big(ids_idx, high_freqs, low_freqs) = b_anti_norm(2, 1);
            end
        end

        toc;
    end
    save([save_dir 'HipopR_cortex_ASB_norm.mat'],"big")
    for ids_idx = 1:length(cortex_id)
        ids = cortex_id(ids_idx);
        fig=figure('Visible','off');

        imagesc(squeeze(big(ids_idx,:,:)));
        colormap hot
        colorbar;
        name=labels{ids};
        title(name);
        clim([0 1])
        saveas(fig,[save_dir name '.png'])
        close all
    end

    % group the ids of regions_cortex that are in the same region
    % and compute the mean of the bicoherence values
    % for each pair of regions
    regions = unique(regions_cortex);

    for region_idx = 1:length(regions)
        region = regions{region_idx};
        ids = find(strcmp(regions_cortex, region));

        % Compute the mean of the bicoherence values
        mean_bicoherence = mean(big(ids, :, :), 1);

        fig=figure('Visible','off');
        imagesc(squeeze(mean_bicoherence));
        colormap hot
        colorbar;
        title(region);
        clim([0 1])
        saveas(fig,[save_dir region '.png'])
        close all
    end

end