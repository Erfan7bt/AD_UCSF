meg_files=dir('../Data/Original/*_meg_rest_60sec.mat');

for i=1:length(meg_files)   
    meg=load([meg_files(i).folder '/' meg_files(i).name]);
    meg=rmfield(meg,'elec');
    meg.hdr=rmfield(meg.hdr, 'elec');
    save(['../Data/Corrected/' meg_files(i).name], 'meg')
end