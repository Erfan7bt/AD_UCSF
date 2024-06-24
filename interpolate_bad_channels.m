function [data_out, bad_channels] = interpolate_bad_channels(rdir, channels, data, channel_layout, channel_name)
% get variance of channels
var_chans = var(data.trial{1,1}(channels,:),0,2);
outliers = isoutlier(var_chans); %three scaled median absolute deviations
if any(outliers)
    bad_channels = channels(outliers);
    cfg = [];
    cfg.channel = channel_name;
    cfg.layout = channel_layout;
    cfg.method = 'triangulation';
    neighbours = ft_prepare_neighbours(cfg, data);
    cfg = [];
    cfg.method = 'weighted';
    % FT bug - badchannel needs to have a name since it's compared to
    % data.label and not just an index
    badchan_cell = data.label(bad_channels);
    cfg.badchannel = badchan_cell;
    cfg.neighbours = neighbours;
    [data_clean] = ft_channelrepair(cfg, data);
    % get cleaned data for particular channel type
    data_out = data_clean.trial{1,1}(channels,:);
else
    data_clean = data;
    data_out = data_clean.trial{1,1}(channels,:);
    bad_channels = [];
end
bad_channels_fig = figure('visible','off');
%subplot(1,2,1)
hold on
plot(1:length(channels),var_chans,'o')
plot(find(outliers),var_chans(outliers),'or')
xlabel('Channels')
ylabel('Variance')
title('Variance of channels')
%subplot(1,2,2)
%ts = reshape(data.trial{1,1},length(data.label),[],100);
%mean_ts = mean(ts,3);
%plot(mean_ts(channels,:)','k')
%plot(mean_ts(bad_channels,:)','r')
%title('Mean of time series')
disp('bad channels:');
fprintf(1, '%s \n', data.label{bad_channels});
saveas(bad_channels_fig,[rdir '/bad_channels.jpg'])
clearvars varchans nbad
