clear all;close all;
%% I/O
infolder='/group_shares/FAIR_ASD/Projects/UCDavis_APP/mapped/connectomes_diffusion_MRtrixCSD_prob';
insubjectlist='NBSprep/SUBS_164ids.txt';
outfolder='NBSprep';

% what parcellation do you want? 'lausanne' or 'GordonParc_reslice'
% lausanne has scale 33, 60, 125, 250, or 500
parc='lausanne';
scale=33; %this only applies to lausanne

%% GO!!
if strcmp(parc,'lausanne')
    cfile=['connectome_scale' num2str(scale) '.mat'];
elseif strcmp(parc,'GordonParc_reslice')
    cfile='connectome_GordonParc_reslice.mat';
end

%% load the list and loop through each subject
fid=fopen(insubjectlist);
slist=textscan(fid,'%s');
fclose(fid);
slist=slist{1,1};
for s=1:length(slist)
    %% load subject's connection matrix into the 3D array
    sub=slist{s};
    load([infolder '/' sub '/connectivity_matrices/' cfile]);
    subject_array_3D(:,:,s)=sc.number_of_fibers;
end

%% save raw connection matrices
savefile=[outfolder '/' num2str(length(slist)) cfile];
save(savefile,'subject_array_3D');

%% only include connections >0 for at least 90% of subjects
non0=sum((subject_array_3D>0),3)>=size(subject_array_3D,3)*.9;
for s=1:length(slist)
    subject_array_3D(:,:,s)=subject_array_3D(:,:,s).*non0;
end
savefile=[outfolder '/exclude0s_' num2str(length(slist)) cfile];
save(savefile,'subject_array_3D');