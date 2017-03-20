clear all;close all;
%% I/O
inmatpref='COMmats/41connectome_'; %COMmats/123connectome_
indesignfile='NBSprep/designANCOVA164_41female.txt';

%subgroups file should have 1's, 2's, 3's, etc. for each subgroup to include, 0 for subjects to exclude
%subgroupfile='explore/traj3groups_maleASD_123connectomesorder_time2only.csv';
%outdesignfile='NBSprep/time2age/traj3groups_maleASD_designANCOVA.txt';

subgroupfile='NBSprep/design_files_behavior/theASD_subgroup_41females.txt';
outdesignfile='NBSprep/design_files_behavior/femaleASD_designANCOVA.txt';

% subgroupfile='explore/lcagroups_maleASD_123connectomesorder.csv';
% outdesignfile='NBSprep/subgroups/lcagroups_maleASD_designANCOVA.txt';

outmatpref='NBSprep/subgroups/femaleASD_COM_ranknormal';

% what parcellation do you want? 'lausanne' or 'GordonParc_reslice'
% lausanne has scale 33, 60, 125, 250, or 500
parc='lausanne';
scale=125; %this only applies to lausanne

% which column of the design matrix to remove?
descolrem = 2; %2 is original diag status; set to 0 not to remove a column

% apply transform to truncated dataset
logtransform=       0; % 0 or 1
ranknormalize=      1; % if 1, overrules logtransform

% exclude missing data? missing entries should be coded -999 in the design file
exclude_missing_data = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GO!!
%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(parc,'lausanne')
    cfile=[inmatpref 'scale' num2str(scale) '.mat'];
    out_cfile=[outmatpref 'scale' num2str(scale) '.mat'];
elseif strcmp(parc,'GordonParc_reslice')
    cfile=[inmatpref 'GordonParc_reslice.mat'];
    out_cfile=[outmatpref 'GordonParc_reslice.mat'];
end

%% load data
load(cfile);
[sx,sy,sz]=size(subject_array_3D);
design=dlmread(indesignfile);
subgroups=dlmread(subgroupfile);
sgs=unique(subgroups(subgroups>0))';
new_array_3D=[];new_design=[];

if length(sgs) > 2
    %% if greater than 2 subgroups, include n-1 extra columns in design matrix, where n = (# of subgroups)
    for g=sgs(1:end-1)
        design=cat(2,design,subgroups==g);
    end
elseif length(sgs) == 2
    %% if exactly 2 subgroups, include one extra column in design matrix for group difference
    extracol = (subgroups==sgs(1))*1 + (subgroups==sgs(2))*-1;
    design=cat(2,design,extracol);  
end

%% truncate conmat and designfile to only include subgroups
for g=sgs
    new_array_3D=cat(3,new_array_3D,subject_array_3D(:,:,subgroups==g));
    new_design=cat(1,new_design,design(subgroups==g,:));
end

%% remove a column from design matrix
if descolrem > 0;
    new_design(:,descolrem)=[];
end

%% find and remove missing data (missing entries should be coded as -999)
[rows_missing,~]=find(new_design==-999);
if length(rows_missing) >= 1
    disp('THERE ARE MISSING DATA HERE!!!');
    disp('THIS WILL OVERWRITE THE OUTPUT DESIGN FILE AND .MAT FILE EXCLUDING THE SUBJECTS WITH MISSING DATA');
    if exclude_missing_data == 0
        disp('EXITING WITH ERROR');
        error('EXITING WITH ERROR. Set exclude_missing_data to 1 and rename the output files');
    end
    new_design(rows_missing,:)=[];
    new_array_3D(:,:,rows_missing)=[];
end

%% rename new array
subject_array_3D=new_array_3D;
%data transformations
subject_array_3D = transform_data(subject_array_3D,logtransform,ranknormalize);

%% save raw connection matrices
save(out_cfile,'subject_array_3D');
dlmwrite(outdesignfile,new_design,' ');