addpath ../fieldtrip-20240704/
ft_defaults
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
segmentedmri = load ('RSID0063_mri_anon.mat');
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
print -dpng natmeg_dip_segorg.png

% use IMDILATE to inflate the segmentation
binary_scalp = imdilate(binary_scalp, strel_bol(1));

figure(2)
tmp = binary_scalp + binary_skull + binary_brain;
imagesc(squeeze(tmp(:,:,128)));

print -dpng natmeg_dip_segdil1.png

% use IMDILATE to inflate the segmentation a bit more
binary_scalp = imdilate(binary_scalp, strel_bol(1));

figure(3)
tmp = binary_scalp + binary_skull + binary_brain;
imagesc(squeeze(tmp(:,:,128)));
print -dpng natmeg_dip_segdil2.png

% revert to the oriiginal binary_scalp
binary_scalp = segmentedmri.scalp + binary_skull;

% use boolean logic together with IMERODE
binary_skull = binary_skull & imerode(binary_scalp, strel_bol(2)); % fully contained inside eroded scalp
binary_brain = binary_brain & imerode(binary_skull, strel_bol(2)); % fully contained inside eroded skull

figure(4)
tmp = binary_scalp + binary_skull + binary_brain;
imagesc(squeeze(tmp(:,:,128)));
print -dpng natmeg_dip_segbool.png
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
figure
ft_plot_mesh(mesh(1), 'facecolor', 'none'); % brain
view([0 -1 0]); % from the right side

figure
ft_plot_mesh(mesh(2), 'facecolor', 'none'); % skull
view([0 -1 0]); % from the right side

figure
ft_plot_mesh(mesh(3), 'facecolor', 'none'); % scalp
view([0 -1 0]); % from the right side

%%
figure
ft_plot_mesh(mesh(1), 'facecolor','r', 'facealpha', 0.4, 'edgecolor', 'none', 'edgealpha', 1);
hold on
ft_plot_mesh(mesh(2), 'facecolor','g', 'facealpha', 0.4, 'edgecolor', 'none', 'edgealpha', 1);
hold on
ft_plot_mesh(mesh(3), 'facecolor','b', 'facealpha', 0.4, 'edgecolor', 'none', 'edgealpha', 1);
%%
% Create a volume conduction model
cfg        = [];
cfg.grad= grad;
cfg.method = 'localspheres'; % You can also specify 'openmeeg', 'bemcp', or another method
cfg.tissue = {'brain'};
headmodel  = ft_prepare_headmodel(cfg, mesh(1));
%%
meg=load('RSID0063_meg_rest_60sec.mat');
grad=meg.grad;
%%
figure
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
%%
cfg=[];
cfg.headmodel=headmodel;
cfg.grad=grad;
sourcemodel=ft_prepare_sourcemodel(cfg);
%%
cfg=[];
cfg.sourcemodel=sourcemodel;
cfg.headmodel= headmodel;
% cfg.normalize = 'yes';
cfg.grad=grad;
lfm=ft_prepare_leadfield(cfg,meg);
