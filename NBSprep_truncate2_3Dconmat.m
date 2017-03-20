clear all;close all;

%% this script takes as input four files: 1) a connection matrix, 2) ANCOVA design file,
%% 3) subject list corresponding to the previous 2 files, and 4) subject list
%% that is a subset of the other subject list. Using just the 2 subject lists... the connection matrix,
%% design file, and original subject list are truncated to the subset list. The truncated subject list may
%% have a different ordering than the input subset list, hence it is outputted as well.

%% I/O
inmatpref='COMmats/164connectome_'; %NBSprep/subgroups/all_maleASD_COM_
indesignANCOVAfile='NBSprep/design_files_behavior/all_maleASD_designANCOVA_beaf_categories.txt';
insubfile='NBSprep/subgroups/SUBS_164_89maleASDids.txt';

%subgroupfile is list of subject IDs that should be a subset of the insubs
subgroupfile='NBSprep/time2age/unordered/SUBS_103_62maleASDids.txt'; %SUBS_103_81maleids.txt, SUBS_103_72ASDids.txt, SUBS_92_74maleids.txt, SUBS_92_64ASDids.txt
outmatpref='NBSprep/time2age/62maleASD_COM'; %81male_COM, 72ASD_COM, 74male_COM, 64ASD_COM

output_new_design_file = 1;
outdesignANCOVAfile='NBSprep/time2age/designANCOVA103_62maleASD_beaf_categories.txt'; %designANCOVA103_81male.txt, designANCOVA103_72ASD.txt, designANCOVA92_74male.txt, designANCOVA92_64ASD.txt
output_new_sub_file = 0;
outsubfile='NBSprep/time2age/SUBS_103_62maleASDids.txt'; %SUBS_103_81maleids.txt, SUBS_103_72ASDids.txt, SUBS_92_74maleids.txt, SUBS_92_64ASDids.txt

% what parcellation do you want? 'lausanne' or 'GordonParc_reslice'
% lausanne has scale 33, 60, 125, 250, or 500
parc='lausanne';
scale=125; %this only applies to lausanne

% apply transform to truncated dataset
logtransform=       0; % 0 or 1
ranknormalize=      1; % if 1, overrules logtransform

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GO!!
%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist(outdesignANCOVAfile, 'file') == 2 && output_new_design_file == 1
    error(['Your outdesignANCOVAfile ' outdesignANCOVAfile ' already exists! Try turning off output_new_sub_file']);
end
if exist(outsubfile, 'file') == 2 && output_new_sub_file == 1
    error(['Your outsubfile ' outsubfile ' already exists! Try turning off output_new_design_file']);
end

if strcmp(parc,'lausanne')
    cfile=[inmatpref 'scale' num2str(scale) '.mat'];
    out_cfile=[outmatpref 'scale' num2str(scale)];
elseif strcmp(parc,'GordonParc_reslice')
    cfile=[inmatpref 'GordonParc_reslice.mat'];
    out_cfile=[outmatpref 'GordonParc_reslice'];
end

%% load supergroup data
load(cfile);
[sx,sy,sz]=size(subject_array_3D);
designANCOVA=dlmread(indesignANCOVAfile);
fileID = fopen(insubfile,'r');
supergroupsubs = textscan(fileID,'%s');
supergroupsubs=supergroupsubs{1};
fclose(fileID);

%% load subgroup data
fileID = fopen(subgroupfile,'r');
subgroupsubs = textscan(fileID,'%s');
subgroupsubs=subgroupsubs{1};
fclose(fileID);

% cycle through subgroupsubs and generate a binary vector of size supergroupsubs
% that indicates which subjects to retain
retainvec=zeros(length(supergroupsubs),1);
for i=1:length(subgroupsubs)
    retainvec(strcmpi(subgroupsubs{i},supergroupsubs))=1;
end

%% truncate conmat, designfile, and sublists to only include subgroup
new_array_3D=subject_array_3D(:,:,retainvec == 1);
new_designANCOVA=designANCOVA(retainvec == 1,:);
new_subs=supergroupsubs(retainvec == 1,:);

%% rename and save new truncated 3-d matrices
subject_array_3D=new_array_3D;
save(out_cfile,'subject_array_3D');
%data transformations
subject_array_3D = transform_data(subject_array_3D,logtransform,ranknormalize);
if logtransform==1 && ranknormalize==0
    save([out_cfile '_logged'],'subject_array_3D');
elseif ranknormalize==1
    save([out_cfile '_ranknormal'],'subject_array_3D');
end

%% write truncated designfile, and sublists
if output_new_design_file == 1
    dlmwrite(outdesignANCOVAfile,new_designANCOVA,' ');
end
if output_new_sub_file == 1
    fid=fopen(outsubfile,'w');
    [nrows,~] = size(new_subs);
    for row = 1:nrows
        fprintf(fid,'%s\n',new_subs{row,:});
    end
    fclose(fid);
end