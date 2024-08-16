%%
% Load the MRI data
addpath ../fieldtrip-20240704/
ft_defaults
mri_files = dir('../Data/Original/*.nii');
%% correct the transfoem

for k = 1:length(mri_files)
mri = ft_read_mri([ mri_files(k).folder  '/'  mri_files(k).name]);
disp('Original Transformation Matrix:');
disp(mri.transform);

mri_corrected = mri;
mri_corrected.transform = [0 -1 0  128.5; 1 0 0 -98.5;0 0 1 -88.5;0 0 0 1];

disp('Corrected Transformation Matrix:');
disp(mri_corrected.transform);

% ft_sourceplot([], mri_corrected);
% ft_sourceplot([], mri);

ft_write_mri(['../Data/Corrected/' mri_files(k).name], mri_corrected,'dataformat','nifti');
end