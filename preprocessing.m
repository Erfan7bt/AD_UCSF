clc
clear all
addpath ../fieldtrip-20230118/
ft_defaults


report_dir = './Results/';
preprocessed_dir = report_dir;
save_prep = preprocessed_dir;
subj='RSID0071';

%% Read and plot raw data 
cfg=[];
data_file= 'RSID0071_meg_rest_60sec.mat';
data=load(data_file);
raw_data=ft_preprocessing(cfg,data);
nchans = length(raw_data.label);

names = 'Axial_GradioMeters';
channel_names=raw_data.label;

cfg=[];
cfg.layout='CTF275_helmet.mat';
layout=ft_prepare_layout(cfg);

plot_data(report_dir, raw_data, layout, channel_names, ['raw_after_import_' names]);

cfg = [];
cfg.resamplefs = 200;
cfg.method = 'downsample';
data_resampled=ft_resampledata(cfg,raw_data);

plot_data(report_dir,data_resampled,layout,channel_names,['resampled_' names])

channel_types = find(~cellfun(@isempty,regexp(data_resampled.label,'^M')));

chanlist = 1:length(data_resampled.label);
all_bad_channels = [];
data_cleaned = data_resampled; % fieldtrip structure

[data_interpolated, bc] = interpolate_bad_channels(report_dir, channel_types, data_resampled, layout, channel_names);
% put into FieldTrip struct
data_cleaned.trial{1,1}(channel_types,:) = data_interpolated;
all_bad_channels = [all_bad_channels;bc];

    %% segment the data into 2s
cfg = [];
cfg.length = 2;
cfg.overlap = 0;
data_seg = ft_redefinetrial(cfg,data_cleaned);
%% reject outlier trials
all_bad_trials = [];
for ch_type = 1:size(channel_types,2)
    bt = detect_bad_trials(report_dir, channel_types(:,ch_type), data_seg);
    all_bad_trials = [all_bad_trials;bt];
end
all_bad_tri = unique(all_bad_trials);
disp('bad trials:');
fprintf(1, '%d \n', all_bad_trials);
% reject trials here and select only meg channels for saving
cfg = [];
cfg.trials = 1:length(data_seg.trial);
cfg.trials(all_bad_tri) = [];
cfg.channel = 'meg';
data_seg = ft_selectdata(cfg,data_seg);
%% plot again for report
for c= 1:1
    plot_data(report_dir, data_seg, layout, channel_names, ['after_rej_' names]);
end
%% save

disp(['This data is saved in ' [save_prep subj] ' under the name' subj]);
ft_write_data([save_prep subj], data_seg, 'dataformat', 'matlab');

