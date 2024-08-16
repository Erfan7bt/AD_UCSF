% bad subjects
%% source rec - subjects:603
% subject 236 (321087) errors during source rec
%% headmodeling - subjects: 623
% checked modeling results
setup
load('/home/space/uniml/datasabzi/results/CamCAN_new/figures/eval_result/eval_result.mat')
no_headmodel = find(any(cellfun(@isempty,squeeze(struct2cell(results)))));
% indices out of 623
% getting bad headmodel subjects out of headmodel check results
bad_hm = [];
for i = 1:length(results)
    a = struct2cell(results(i));
    b = [a{:}];
    if ~isempty(find(~b,1))
        bad_hm = [bad_hm, i];
    end
end
%%
% manually compiled from headmodel check
bad_headmodel = [23, 65, 78, 115, 123, 255, 282, 309, 322, 345, 352, 385, 451, 468, 474, 535, 589, 592]; 
bad_subjects_hm = unique([no_headmodel,bad_headmodel]);
% find subject ids from folder names
hm_dir = dir('/home/space/uniml/datasabzi/results/CamCAN_old/figures/headmodeling/');
hm_dir = hm_dir(3:end);
bad_hm_ids = {hm_dir(bad_subjects_hm).name};
% convert to number
for i=1:length(bad_hm_ids)
    bad_hm_ids{i} = str2double(bad_hm_ids{i}(7:end));
end
bad_hm_ids = [bad_hm_ids{:}];
% some of these aren't in the MEG subjects, remove
in_meg = ismember(bad_hm_ids,sub_info(:,1));
bad_hm_ids(~in_meg) = [];
%% MEG artefacts after prep - subject: 631
% checked topoplots for artefacts
bad_meg = [10, 59, 70, 190, 236, 247, 261, 267, 272, 279, 286, 299, 303, 326, 356, 381, 382, 385, 389, 403, 438, 461, 481, 515, 528, 529, 543, 565, 579, 598, 619, 622];
bad_meg_ids = sub_info(bad_meg,1)';
% combine
bad_subjects_all = unique([bad_meg_ids,bad_hm_ids]);
bad_sub_id = find(ismember(sub_info(:,1),bad_subjects_all));
%%
good_sub = ~ismember(sub_info(:,1),bad_subjects_all);
sub_info_clean = sub_info(good_sub,:);
%save([repo 'sub_info_clean.mat'],'sub_info_clean')