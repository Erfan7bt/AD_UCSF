% load anatomical data from the brainstorm_db folder, convert to MNI,
% extrapolate to high res, make leadfield, save to bs_results for each
% subject
anat_folder = '/home/erfan/Thesis/brainstorm_db/AD_freesurfer/anat/';
data_folder = '/home/erfan/Thesis/brainstorm_db/AD_freesurfer/data/';
AD_dir = '/home/erfan/Thesis/ADanonShare/';

bs_dir = dir([anat_folder 'RSID*']);

bs_names = {bs_dir.name}';

for n = 1:length(bs_names)
    % based on proc_bs_files.m
    disp(['loading data for subject ' bs_names{n}])
    result_folder = [AD_dir '/Results/headmodeling/' bs_names{n} '/'];

    f = dir([anat_folder bs_names{n}]);
    if length(f) >= 20 %check if all files are present
        if ~isfolder(result_folder)
            mkdir(result_folder)
            %check for all files
            Mixed_cortex_lowres = load([anat_folder bs_names{n} '/tess_concat.mat']);
            nvox = size(Mixed_cortex_lowres.Vertices, 1);
            % highres
            Mixed_cortex_highres = load([anat_folder bs_names{n} '/tess_concat_03.mat']);
            % lowres
            Mixed_cortex = load([anat_folder bs_names{n} '/tess_concat_02.mat']);
            % convert to MNI
            disp('converting to MNI')
            Mixed_cortex = conv_mni(Mixed_cortex);
            Mixed_cortex_highres = conv_mni(Mixed_cortex_highres);
            Mixed_cortex_lowres = conv_mni(Mixed_cortex_lowres);
            % calculate extrapolations
            disp('extrapolations to high resolution cortex')
            mi = [];
            in_normal_to_high = [];
            for ii = 1:size(Mixed_cortex_highres.Vertices, 1)
                [mi(ii), in_normal_to_high(ii)] = min(eucl(Mixed_cortex_highres.Vertices(ii, :), Mixed_cortex.Vertices));
            end
            mi = [];
            in_low_to_high = [];
            for ii = 1:size(Mixed_cortex_highres.Vertices, 1)
                [mi(ii), in_low_to_high(ii)] = min(eucl(Mixed_cortex_highres.Vertices(ii, :), Mixed_cortex_lowres.Vertices));
            end
            [~, ia, ib] = intersect(Mixed_cortex.Vertices, Mixed_cortex_lowres.Vertices, 'rows');
            [~, ic] = sort(ib);
            in_normal_to_low = ia(ic);
            % load BEM
            disp('loading BEM model')
            headmodel = load([data_folder bs_names{n} '/@default_study/headmodel_mix_openmeeg.mat']);
            leadfield = permute(reshape(headmodel.Gain, [], 3, length(headmodel.GridLoc)), [1 3 2]);
            disp('saving result')
            save([result_folder 'bs_results'], 'Mixed_cortex', 'Mixed_cortex_highres', 'Mixed_cortex_lowres', 'leadfield', ...
                'in_normal_to_high', 'in_low_to_high', 'in_normal_to_low');
            clearvars ia ib ic ii mi so
        end
     else
         disp(['no data for subject ' num2str(n)])
         continue
    end
end
function cort = conv_mni(cort)
cort.Vertices = cort.Vertices(:, [2 1 3]);
cort.Vertices(:, 1) = -cort.Vertices(:, 1);
end
