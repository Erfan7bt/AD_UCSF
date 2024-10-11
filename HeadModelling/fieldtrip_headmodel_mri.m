
[ret, name] = system('hostname');

if startsWith(name,'ra')

    addpath '/data/erfan/brainstorm3'

    fs_folder = '/data/erfan/freesurfer/subjects/';
    
    BrainstormDbDir = '/data/erfan/brainstorm_db/';
    
    AD_dir = '/home/erfan/AD_UCSF/';
else 

    fs_folder = '/home/erfan/Thesis/freesurfer/subjects/';
    
    BrainstormDbDir = '/home/erfan/Thesis/brainstorm_db/';
    
    AD_dir = '/home/erfan/Thesis/ADanonShare/';
end 

fig_folder = [AD_dir '/Results/Figures/'];

addpath /home/erfan/Thesis/fieldtrip-20240704/

Data_dir=[AD_dir '/Data/Original/'];

ft_defaults
datafiles=dir([Data_dir '*mri_anon.mat'])
%% BEM head model from original mri and segmentation from scratch 
% mri=ft_read_mri('RSID0063_mri_anon.nii');
% 
% cfg = [];
% cfg.method = 'interactive';
% cfg.coordsys = 'ctf'; % the desired coordinate system
% mri_realigned = ft_volumerealign(cfg, mri);
% ft_sourceplot([],mri)
% ft_sourceplot([],mri)
% %%
% % mri_resliced=mri_realigned;
% cfg = [];
% cfg.method = 'spline';
% mri_resliced = ft_volumereslice(cfg, mri_realigned);
% 
% cfg = [];
% cfg.method = 'ortho';
% ft_sourceplot(cfg, mri_resliced)
% ft_determine_coordsys(mri_resliced)
%%
% provided segmentation
% loop over all the files in the directory that end with mri_anon.mat
for i=1:length(datafiles)
    sbj=datafiles(i).name
segmentedmri = load ([Data_dir sbj]);
ft_checkdata(segmentedmri, 'feedback', 'yes') % display some information about the segmentation

% 
% cfg           = [];
% cfg.output    = {'brain', 'skull', 'scalp'};
% % cfg.spmmethod = 'new';
% cfg.brainsmooth = 7;
% cfg.scalpsmooth = 7;
% cfg.skullsmooth = 7;
% segmentedmri  = ft_volumesegment(cfg, mri_resliced);

segmentedmri_indexed = ft_checkdata(segmentedmri, 'segmentationstyle', 'indexed')

% segmentedmri_indexed.anatomy = mri_resliced.anatomy;

cfg = [];
cfg.method = 'ortho';
cfg.anaparameter = 'anatomy';
cfg.funparameter = 'tissue';
cfg.funcolormap = [
  0 0 0
  1 0 0
  0 1 0
  0 0 1
  ];
ft_sourceplot(cfg, segmentedmri_indexed)
% cfg = [];
% cfg.method = 'ortho';
% ft_sourceplot(cfg, mri)
%%
binary_brain = segmentedmri.brain;
binary_skull = segmentedmri.skull | binary_brain;
binary_scalp = segmentedmri.scalp | binary_brain | binary_skull;
%%
close all

% using ft_sourceplot I identified the cross-section with voxel
% indices [107 100 100] where the problem is visible and I will
% plot that intersection multiple times

figure(1)
tmp = binary_scalp + binary_skull + binary_brain;
imagesc(squeeze(tmp(:,:,128)));
% print -dpng natmeg_dip_segorg.png

% use IMDILATE to inflate the segmentation
binary_scalp = imdilate(binary_scalp, strel_bol(1));

figure(2)
tmp = binary_scalp + binary_skull + binary_brain;
imagesc(squeeze(tmp(:,:,128)));

% print -dpng natmeg_dip_segdil1.png

% use IMDILATE to inflate the segmentation a bit more
binary_scalp = imdilate(binary_scalp, strel_bol(1));

figure(3)
tmp = binary_scalp + binary_skull + binary_brain;
imagesc(squeeze(tmp(:,:,128)));
% print -dpng natmeg_dip_segdil2.png

% revert to the oriiginal binary_scalp
binary_scalp = segmentedmri.scalp + binary_skull;

% use boolean logic together with IMERODE
binary_skull = binary_skull & imerode(binary_scalp, strel_bol(2)); % fully contained inside eroded scalp
binary_brain = binary_brain & imerode(binary_skull, strel_bol(2)); % fully contained inside eroded skull

figure(4)
tmp = binary_scalp + binary_skull + binary_brain;
imagesc(squeeze(tmp(:,:,128)));
% print -dpng natmeg_dip_segbool.png
close all
%%
mri_segmented2 = segmentedmri;
% insert the updated binary volumes, taking out the center part for skull and scalp
mri_segmented2.brain    = binary_brain;
mri_segmented2.skull    = binary_skull & ~binary_brain;
mri_segmented2.scalp    = binary_scalp & ~binary_brain & ~binary_skull;
mri_segmented2.combined = binary_scalp + binary_skull + binary_brain; % only for plotting
%%
cfg = [];
cfg.funparameter = 'combined';
cfg.funcolormap = 'jet';
ft_sourceplot(cfg, mri_segmented2);
% this has to be removed, otherwise ft_prepare_mesh gets confused
mri_segmented2 = rmfield(mri_segmented2, 'combined');
%%
cfg = [];
cfg.tissue = {'brain', 'skull', 'scalp'};
cfg.method = 'projectmesh'; 
cfg.numvertices = [3000 900 900];
mesh = ft_prepare_mesh(cfg, mri_segmented2);
% cfg = [];
% cfg.method = 'projectmesh';
% cfg.tissue = 'brain';
% cfg.numvertices = 3000;
% mesh_brain = ft_prepare_mesh(cfg, mri_segmented2);
% 
% cfg = [];
% cfg.method = 'projectmesh';
% cfg.tissue = 'skull';
% cfg.numvertices = 2000;
% mesh_skull = ft_prepare_mesh(cfg, mri_segmented2);
% 
% cfg = [];
% cfg.method = 'projectmesh';
% cfg.tissue = 'scalp';
% cfg.numvertices = 1000;
% mesh_scalp = ft_prepare_mesh(cfg, mri_segmented2);
%%
% figure
% ft_plot_mesh(mesh(1), 'facecolor', 'none'); % brain
% view([0 -1 0]); % from the right side
% 
% figure
% ft_plot_mesh(mesh(2), 'facecolor', 'none'); % skull
% view([0 -1 0]); % from the right side
% 
% figure
% ft_plot_mesh(mesh(3), 'facecolor', 'none'); % scalp
% view([0 -1 0]); % from the right side

%%
% figure
% ft_plot_mesh(mesh(1), 'facecolor','r', 'facealpha', 0.4, 'edgecolor', 'none', 'edgealpha', 1);
% hold on
% ft_plot_mesh(mesh(2), 'facecolor','g', 'facealpha', 0.4, 'edgecolor', 'none', 'edgealpha', 1);
% hold on
% ft_plot_mesh(mesh(3), 'facecolor','b', 'facealpha', 0.4, 'edgecolor', 'none', 'edgealpha', 1);
%%
sbj=strrep(sbj,'_mri_anon.mat','');
meg=load([Data_dir sbj '_meg_rest_60sec.mat']);
grad=meg.grad;
%%
sensor=figure('visible','off');
ft_plot_sens(grad,'axis',true, 'unit', 'cm')
hold on 
% ft_plot_headmodel(headmodel)
ft_plot_mesh(mesh(1), 'facecolor','r', 'facealpha', 0.4, 'edgecolor', 'none', 'edgealpha', 1);
hold on
ft_plot_mesh(mesh(2), 'facecolor','g', 'facealpha', 0.4, 'edgecolor', 'none', 'edgealpha', 0.1);
hold on
ft_plot_mesh(mesh(3), 'facecolor','b', 'facealpha', 0.4, 'edgecolor', 'none', 'edgealpha', 0.1);
% ft_plot_headmodel(headmodel, 'facealpha', 0.5, 'edgecolor', 'none'); % "lighting phong" does not work with opacity
material dull
camlight
saveas(sensor, [fig_folder sbj '/ft_sensors_on_orig_head_1.png']);
az = 0; el = 0; view(az, el)
saveas(sensor, [fig_folder sbj '/ft_sensors_on_orig_head_2.png']);

az = -90; el = 0; view(az, el)
saveas(sensor, [fig_folder sbj '/ft_sensors_on_orig_head_3.png']);

az = -180; el = 0; view(az, el)
saveas(sensor, [fig_folder sbj '/ft_sensors_on_orig_head_4.png']);

az = -270; el = 0; view(az, el)
saveas(sensor, [fig_folder sbj '/sensors_ft_on_orig_head_5.png']);

close all
end