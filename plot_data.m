function plot_data(sub_report, data, layout, channels, desc)
if length(data.trial) == 1
    % segment the data into 2s
    cfg = [];
    cfg.length = 2;
    cfg.overlap = 0;
    data_epoched = ft_redefinetrial(cfg,data);
else
    data_epoched = data;
end
%run FFT on the data
cfg = [];
cfg.channel = channels;
cfg.method = 'mtmfft';
cfg.output = 'pow';
cfg.taper = 'dpss';
cfg.pad = 'nextpow2'; % improves speed
cfg.tapsmofrq = 2;
cfg.foilim = [0.5 55];
% avg of all trials, PSD
freq_seg=ft_freqanalysis(cfg,data_epoched);
cfg.layout = layout;
psd = figure('visible','off');
semilogy(freq_seg.freq,freq_seg.powspctrm)
grid on
xlim([0.5 55]);
xlabel('frequency')
ylabel('power')
topo = figure('visible','off');
% find the indices of frequencies for the delta, theta, alpha and beta
% bands
freq_bounds = {find(freq_seg.freq<4); ...
    find(freq_seg.freq>4 & freq_seg.freq<7); ...
    find(freq_seg.freq>8 & freq_seg.freq<13); ...
    find(freq_seg.freq>14 & freq_seg.freq<30)};
fnames = {'delta', 'theta', 'alpha', 'beta'};
%% plot them in subplots
cfg = [];
cfg.figure = 'gcf';
cfg.layout = layout;
cfg.channel = channels;
for nfreq = 1:length(freq_bounds)
    temp_freq = freq_seg;
    subplot(2,2,nfreq)
    temp_freq.freq = temp_freq.freq(freq_bounds{nfreq});
    temp_freq.powspctrm = temp_freq.powspctrm(:,freq_bounds{nfreq});
    ft_topoplotER(cfg, temp_freq); colorbar
    title(fnames{nfreq})
end
if ~isfolder(sub_report)
    mkdir(sub_report)
end
%%
saveas(topo,[sub_report '/topo_' desc '.png'])
saveas(psd,[sub_report '/spectrum_' desc '.png'])
end