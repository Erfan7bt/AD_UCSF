
settings.target_fs = 100; %sampling frequency
settings.lcmv_reg = 0.05;
settings.fres = settings.target_fs;
settings.smoothcortex = 0.3; %less smoothing
% epoch length
settings.winlength = 2*settings.target_fs; 
% model order for the spectral/connectivity analysis
settings.morder = 20;
load cm17
settings.cma = cm17a;
% how many PCA components to keep for each ROI
settings.nPCA = 3; 

%% loop over subjects
nsub = 1
% load MEG data
tic
subj = 'sub-CC110037';
disp(['processing subject ' subj])
load([subj '.mat']);
nchan = size(dat.label,1);
ntri = size(dat.trial,2); 
nsam = size(dat.trial{1,1},2);
% load Brainstorm files

try
    load(['./bs_results.mat']);
    disp('loading data')
catch
    % skip if no data
    disp('no Brainstorm data found')
end
% put the data in the right format (nchan x nsam x ntri):
data = cat(3,dat.trial{1,:});
% upscale for numerical stability
    data = data.*10^10;
    leadfield = leadfield.*10^10;
    leadfield = leadfield(1:nchan, :, :);
% change the leadfield to remove radial dimension for MEG using the function from Franziska
leadfield = fp_get_lf(leadfield);
% frequencies with 0.5 Hz resolution
frqs = linspace(0, settings.target_fs/2, nsam/2+1)';
% number of bootstrap samples
nbootstrap = 0;
% number of ROIs in the Desikan-Kiliany Atlas
nROI = length(cortex.Atlas(3).Scouts);
% ROI labels
labels = {cortex.Atlas(3).Scouts.Label};
ind_cortex = [];
ind_roi = {};
for iROI = 1:nROI
  ind_roi{iROI} = cortex.Atlas(3).Scouts(iROI).Vertices;
  ind_cortex = cat(1, ind_cortex, ind_roi{iROI});
  [~, ind_roi_cortex{iROI}, ~] = intersect(ind_cortex, ind_roi{iROI});
end
nvox = length(ind_cortex);
leadfield = leadfield(:, ind_cortex, :);
%%
% get cross-spectrumtt
disp('computing cross-spectrum')
conn = data2spwctrgc(data, settings.fres, settings.morder, 0, nbootstrap, [], {'CS'});
CS_sensor = conn.CS; 
disp('applying LCMV')
% LCMV projection kernel
C = cov(data(:,:)'); 
alpha = settings.lcmv_reg*trace(C)/length(C);
Cr = C + alpha*eye(nchan);
[~, P] = lcmv_meg(Cr, leadfield, struct('alpha', 0, 'onedim', 0));
% project to source space
% 2 dimensions for MEG
source_voxel_data = reshape(data(:, :)'*P(:, :), nsam*ntri, nvox, 2); 
source_voxel_data = 10^3*source_voxel_data;
% keep only the first nPCA strongest components for each ROI
% initialize
disp('computing source power')
source_roi_data = []; 
source_roi_power = [];
source_roi_power_norm = [];
source_roi_power_total = [];
source_roi_power_total_norm = [];
varex = ones(nROI, size(dat.label,1));
nPCAs = [];
%%
for iROI = 1:nROI
  data_ = source_voxel_data(:, ind_roi_cortex{iROI}, :);
  va_ = var(data_(:, :));
  source_roi_power(iROI) = sum(va_)';
  source_roi_power_norm(iROI) = source_roi_power(iROI)/length(ind_roi_cortex{iROI});
  % updated computation
  P_ = P(:, ind_roi_cortex{iROI}, :);
  P_ = P_(:, :);
  CS_ROI = [];
  for ifreq = 1:size(CS_sensor, 1)
    CS_ROI(ifreq, :, :) = P_'*squeeze(CS_sensor(ifreq, :, :))*P_;
  end
  source_power_all{iROI} = abs(cs2psd(CS_ROI));
  source_roi_power_total(:, iROI) = sum(abs(cs2psd(CS_ROI)), 2); 
  source_roi_power_total_norm(:, iROI) = source_roi_power_total(:, iROI)/length(ind_roi_cortex{iROI}); 
  % optional z-scoring - removed
  % data_(:, :) = zscore(data_(:, :))*sqrt(mean(va_));
  % PCA/SVD
  [data_, S_] = svd(data_(:, :), 'econ');
  % variance explained
  vx_ = cumsum(diag(S_).^2)./sum(diag(S_).^2);
  invx = 1:min(length(vx_), size(dat.label,1));
  varex(iROI, invx) = vx_(invx);
  if settings.nPCA < 1
    nPCAs(iROI) = min(find(varex(iROI, invx) > settings.nPCA)); 
  else
    nPCAs(iROI) = settings.nPCA;
  end
  % keep nPCA components with correct unit and scaling
  source_roi_data = cat(2, source_roi_data, data_(:, 1:nPCAs(iROI))*S_(1:nPCAs(iROI), 1:nPCAs(iROI)));
end
%% get data and indices for saving

% permute data
source_roi_data = permute(reshape(source_roi_data, nsam, ntri, []), [3 1 2]);

% indices for PCAs
beg_inds = cumsum([1 nPCAs(1:end-1)]);
end_inds = cumsum([nPCAs]);
% ranges
for iroi = 1:nROI
  PCA_inds{iroi} = beg_inds(iroi):end_inds(iroi);
end
% combinations of nPCAs to use in connectivity computations (e.g. [1 2 3]
% -> [4 5 6]
inds = {}; ninds = 0;
for iroi = 1:nROI
  for jroi = (iroi+1):nROI
    inds{ninds+1} = {(iroi-1)*settings.nPCA + [1:settings.nPCA], (jroi-1)*settings.nPCA + [1:settings.nPCA]};    
    %inds{ninds+2} = {(jroi-1)*settings.nPCA + [1:settings.nPCA], (iroi-1)*settings.nPCA + [1:settings.nPCA]};  
    ninds = ninds + 1; 
  end
end

if settings.nPCA < 1
    inds = PCA_inds; 
end

% result folder
result_folder_sub = ['./' subj];

% save results 
disp('saving results')
save([result_folder_sub '/source_rec_results.mat'], 'source_roi_power', 'source_roi_power_norm', ...
    'source_roi_power_total','source_roi_power_total_norm','source_power_all', 'conn', 'settings','source_roi_data','inds','varex', ...
     'ind_cortex', 'nchan', 'nPCAs', 'beg_inds', 'end_inds', 'PCA_inds','nbootstrap');
toc
