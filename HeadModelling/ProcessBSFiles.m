% load anatomical data from the brainstorm_db folder, convert to MNI,
% extrapolate to high res, make leadfield, save to bs_results for each
% subject
anat_folder = '/datasabzi/results/CamCAN_old/brainstorm_db/CamCAN/anat/';
data_folder = '/datasabzi/results/CamCAN_old/brainstorm_db/CamCAN/data/';
bs_dir = dir([anat_folder '*CC*']);
bs_names = {bs_dir.name}';
for n = 1:length(bs_names)
    % based on proc_bs_files.m
    disp(['loading data for subject ' bs_names{n}])
    result_folder = ['/datasabzi/results/CamCAN_new/bs_processed/' bs_names{n} '/'];
    f = dir([anat_folder bs_names{n}]);
    if length(f) >= 20 %check if all files are present
        if ~isfolder(result_folder)
            mkdir(result_folder)
            %check for all files
            cortex = load([anat_folder bs_names{n} '/tess_cortex_mid_low.mat']);
            nvox = size(cortex.Vertices, 1);
            % highres
            cortex_highres = load([anat_folder bs_names{n} '/tess_cortex_mid_high.mat']);
            % lowres
            cortex_lowres = load([anat_folder bs_names{n} '/tess_cortex_mid_low_2000V.mat']);
            % convert to MNI
            disp('converting to MNI')
            cortex = conv_mni(cortex);
            cortex_highres = conv_mni(cortex_highres);
            cortex_lowres = conv_mni(cortex_lowres);
            % calculate extrapolations
            disp('extrapolations to high resolution cortex')
            mi = [];
            in_normal_to_high = [];
            for ii = 1:size(cortex_highres.Vertices, 1)
                [mi(ii), in_normal_to_high(ii)] = min(eucl(cortex_highres.Vertices(ii, :), cortex.Vertices));
            end
            mi = [];
            in_low_to_high = [];
            for ii = 1:size(cortex_highres.Vertices, 1)
                [mi(ii), in_low_to_high(ii)] = min(eucl(cortex_highres.Vertices(ii, :), cortex_lowres.Vertices));
            end
            [~, ia, ib] = intersect(cortex.Vertices, cortex_lowres.Vertices, 'rows');
            [~, ic] = sort(ib);
            in_normal_to_low = ia(ic);
            % load BEM
            disp('loading BEM model')
            headmodel = load([data_folder bs_names{n} '/@default_study/headmodel_surf_openmeeg.mat']);
            leadfield = permute(reshape(headmodel.Gain, [], 3, nvox), [1 3 2]);
            disp('saving result')
            save([result_folder 'bs_results'], 'cortex', 'cortex_highres', 'cortex_lowres', 'leadfield', ...
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
