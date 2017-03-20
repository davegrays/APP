clear;clc;close all;
addpath(genpath('/home/david/Projects/matlab_scripts'));

%% USER DEFINED INPUTS
ANCOVAdesign=       load('NBSprep/designANCOVA164_123male.txt');
groups=             load('explore/groups_164ids_male_diagnoses.txt');
inmat=              'NBSprep/123connectome_scale125.mat';
ROI_info_file=      'NBSprep/parcellation_info/ROI_INFO_scale125.txt';

logtransform=       0; % 0 or 1
ranknormalize=      1; % if 1, overrules logtransform

remove_lastnode=    1; %0 or 1 to get rid of brainstem

regress_nuisance=   1;
% set nuisance covariates here
age=ANCOVAdesign(:,3);
TCV=ANCOVAdesign(:,4);

% set effect of interest here (diagnosis, behavior, etc.)
effect=ANCOVAdesign(:,2);

%% BEGIN SCRIPT
%% load data
load(inmat);
[sx, sy, sz]=size(subject_array_3D);
orig_array_3D=subject_array_3D;
nuisance_predictors=[ones(sz,1) age TCV];
all_predictors=[nuisance_predictors effect];
fid=fopen(ROI_info_file);
ROI_info=textscan(fid,'%d %d%d%d %s %s%s %s%s %s%s','emptyValue',0,'HeaderLines',1,'TreatAsEmpty','na');
fclose(fid);
labels=ROI_info{1,5};
    
%% get rid of last node (brainstem) so that matrix has hemispheric symmetry 
if remove_lastnode == 1
    labels(end)=[];
    subject_array_3D(sx,:,:)=[];
    subject_array_3D(:,sx,:)=[];
    [sx, sy, sz]=size(subject_array_3D); %reset matrix size
end

% %% extract ipsis, homos, heteros, etc.
% %node ordering of input matrix is supposed to be N regions on the left followed by N regions on the right
% %(Lausanne parcellations are actually the opposite, so I switched the naming of the outputs here)
%[M_r,M_l,I,hom,contras_r,contras_l] = ipsi(subject_array_3D);

[Gout G_METRICS]=graphtheory_loc_3D(subject_array_3D,'wei',1,0,0,0,0,0,0,0,0,0,0,0,0,0); %STREN
%[Gout G_METRICS]=graphtheory_loc_3D(subject_array_3D,'wei',0,0,0,0,0,0,1,0,0,0,0,0,0,0); %EIG CEN
%[Gout G_METRICS]=graphtheory_loc_3D(subject_array_3D,'wei',0,0,0,0,0,0,0,1,0,0,0,0,0,0); %TOTAL COM
%[Gout G_METRICS]=graphtheory_loc_3D(subject_array_3D,'wei',0,0,0,0,0,0,0,0,1,0,0,0,0,0); %SUBGRAPH CEN
%[Gout G_METRICS]=graphtheory_loc_3D(subject_array_3D,'wei',0,0,0,0,0,0,0,0,0,1,0,0,0,0); %CBC

%[Gout G_METRICS]=graphtheory_glob_3D(subject_array_3D,'wei',0,0,0,0,0,1,0,0); %

%% set your outcome variable inside the squeeze function (outcome var is matrix of size # subjects X # nodes)
%hom, contras_r, contras_l, or Gout for the graph metric
brain_nodes = squeeze(Gout);

%% loop through nodes and compute ANCOVA t-stats
for i=1:size(brain_nodes,2)
    outcome = brain_nodes(:,i);
    
    if ranknormalize == 1
        rank = tiedrank(outcome);
        p = rank / ( length(rank) + 1 ); %# +1 to avoid Inf for the max point
        outcome = norminv( p, 0, 1 ) * std(outcome) + mean(outcome);
    end
    
    %% regress out nuisance effects
    % b=nuisance_predictors\outcome;
    % regressed_y=outcome-nuisance_predictors*b + mean(outcome);
    % %% test effect of interest
    % [H P C T_STATS] = ttest2(regressed_y(effect==1), regressed_y(effect==-1),'vartype','unequal');
    % P
    % T_STATS

    if regress_nuisance == 1
        fit=fitlm([age TCV effect],outcome,'CategoricalVars',3); %essentially an ANCOVA if you take the tstat of the final predictor
        nodes(i,:)=round(10000*[fit.Coefficients.tStat(4) fit.Coefficients.pValue(4)])/10000;
    else
        fit=fitlm(effect,outcome,'CategoricalVars',1);
        nodes(i,:)=round(10000*[fit.Coefficients.tStat(2) fit.Coefficients.pValue(2)])/10000;
    end
    tstats(i)=nodes(i,1);
    pvals(i)=nodes(i,2);
end

FDR = mafdr(pvals,'BHFDR','True');
disp(['FDR-corrected P = ' num2str(min(FDR))]);

siglabels=labels(pvals<0.05);
signodes=nodes(pvals<0.05,:);

fprintf('\nASD greater than TYP\n')
for n=1:length(siglabels)
    if signodes(n,1) < 0
        disp([siglabels{n} ', tstat = ' num2str(signodes(n,1), '%0.3f') ', p = ' num2str(signodes(n,2), '%0.4f')]);
    end
end
fprintf('\nTYP greater than ASD\n')
for n=1:length(siglabels)
    if signodes(n,1) > 0
        disp([siglabels{n} ', tstat = ' num2str(signodes(n,1), '%0.3f') ', p = ' num2str(signodes(n,2), '%0.4f')]);
    end
end

% nodestrengths=mean(brain_nodes,1);
% [y,i]=sort(nodestrengths);
% labels(i)