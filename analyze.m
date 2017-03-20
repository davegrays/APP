clear;clc;close all;
addpath(genpath('/home/david/Projects/matlab_scripts'));

%% USER DEFINED INPUTS
% ANCOVAdesign=       load('NBSprep/subgroups/lcagroups_maleASD_designANCOVA.txt');
% groups=             load('explore/groups_maleASD_lcagroups.txt');
% inmat=              'NBSprep/subgroups/lcagroups_maleASD_COM_ranknormal_scale125.mat'; 
% ANCOVAdesign=       load('NBSprep/designANCOVA164_123male_plusDQ.txt');
% groups=             load('explore/groups_164ids_male_diagnoses.txt');
% inmat=              'COMmats/123connectome_scale125.mat';
% ANCOVAdesign=       load('NBSprep/designANCOVA164_108ASD.txt');
% groups=             load('explore/groups_164ids_ASD_genders.txt');
% inmat=              'COMmats/108connectome_scale125.mat';
% ANCOVAdesign=       load('NBSprep/designANCOVA164_41female.txt');
% groups=             load('explore/groups_164ids_female_diagnoses.txt');
% inmat=              'COMmats/41connectome_scale125.mat';
ANCOVAdesign=       load('NBSprep/designANCOVA164_56TYP.txt');
groups=             load('explore/groups_164ids_TYP_genders.txt');
inmat=              'COMmats/56connectome_scale125.mat';
% ANCOVAdesign=       load('NBSprep/designANCOVA164_164_malethenfemale.txt');
% groups=             load('explore/groups_164ids_malethenfemale.txt');
% inmat=              'COMmats/164connectome_malethenfemale_scale125.mat';
% ANCOVAdesign=       load('NBSprep/design_files_behavior/all_maleASD_designANCOVA_megalen_categories_noTCV.txt');
% groups=             load('explore/groups_maleASD_megalen.txt');
% inmat=              'NBSprep/subgroups/all_maleASD_COM_scale125.mat';

ROIinfofile_int=     'NBSprep/parcellation_info/ROI_INFO_int_scale125_rescaled.txt';

logtransform=       0; % 0 or 1
ranknormalize=      1; % if 1, overrules logtransform

run_mds_svd=        1; %0 or 1
caret_node_stren=   0; %0 or 1
    stren_perc=     0.2;
    outfilename=    'explore/caret_files/scale125_NodeStrengths';
NBS_vs_node_stren=  0; %0 or 1
    nbs_effects=    'NBSresults/caret_files/scale125_COM_t-test_thresh2.6_p0.012.POINTSZ.mat';

% set nuisance covariates here
age=ANCOVAdesign(:,3);
TCV=ANCOVAdesign(:,4);
%DQ=ANCOVAdesign(:,5);

predictors=[ones(length(age),1) age TCV]; %can decide here whether to add TCV or DQ or not

%% BEGIN SCRIPT
%load data
load(inmat);
[sx, sy, sz]=size(subject_array_3D);
orig_array_3D=subject_array_3D;

%data transformations
subject_array_3D = transform_data(subject_array_3D,logtransform,ranknormalize);

%% zero everything below upper triangle
for s=1:sz;subject_array_3D(:,:,s)=triu(subject_array_3D(:,:,s),1);end

%% convert elements in upper triangle to residuals after nuisance regression
B1=zeros(sx,sx);B2=zeros(sx,sx);
disp('nuisance regression');
for i=1:sx
    for j=i+1:sx
        outcome=subject_array_3D(i,j,:);outcome=outcome(:);outcome_mean=mean(outcome);
        [coefs,~,residuals] = regress(outcome,predictors); %nuisance regression
        subject_array_3D(i,j,:)=residuals+outcome_mean; %add mean back to each variable
        B1(i,j)=coefs(1);B2(i,j)=coefs(2); %regression coefficients here for inspection
    end
end

%% re-symmetrize the residuals matrix
for s=1:sz;subject_array_3D(:,:,s)=subject_array_3D(:,:,s)+subject_array_3D(:,:,s)';end

if run_mds_svd==1
    %% produce mds and svd figs on UNREGRESSED intersubject corrmat 
    %% use all edges for mds
    disp('visualizing unregressed data');
    [coords_orig,p_intragroupsimilarity_orig]=matrices2multidim(orig_array_3D,groups,2,'Pearson',0,1);
    hold on;
    for i=unique(groups)'
        groupi=(groups==i);
        ellipse(std(coords_orig(groupi,1)),std(coords_orig(groupi,2)),0,mean(coords_orig(groupi,1)),mean(coords_orig(groupi,2)),'k');
        plot(mean(coords_orig(groupi,1)),mean(coords_orig(groupi,2)),'k+');
    end
    hold off;
    figure();
    %% use node strengths for svd (edges too slow)
    matrices2svdscores(orig_array_3D,groups,1);

    figure();
    %% produce mds and svd figs on REGRESSED intersubject corrmat 
    %% use all edges for mds
    disp('visualizing regressed data');
    [coords_regressed,p_intragroupsimilarity_nuisanceregressed]=matrices2multidim(subject_array_3D,groups,2,'Pearson',0,1);
    hold on;
    for i=unique(groups)'
        groupi=(groups==i);
        ellipse(std(coords_regressed(groupi,1)),std(coords_regressed(groupi,2)),0,mean(coords_regressed(groupi,1)),mean(coords_regressed(groupi,2)),'k');
        plot(mean(coords_regressed(groupi,1)),mean(coords_regressed(groupi,2)),'k+');
    end
    hold off;
    figure();
    %% use node strengths for svd (edges too slow)
    matrices2svdscores(subject_array_3D,groups,1);
end

if caret_node_stren == 1
    %% generate caret brain of mean node strengths or mean node Gtots
    meanmat=mean(subject_array_3D,3);
    meanstrengths=mean(meanmat);
    sv=sort(meanstrengths);
    min_stren=sv(round((1-stren_perc)*length(sv)));
    meanstrengths(meanstrengths<min_stren)=0;
    caretfoci_filemaker(ROIinfofile_int,meanstrengths,ones(length(meanstrengths),1),outfilename,'hold')
end
p_intragroupsimilarity_orig
p_intragroupsimilarity_nuisanceregressed

if NBS_vs_node_stren == 1
    %% look at association between NBS node effects and mean node strengths (or mean node Gtots)
    meanmat=mean(subject_array_3D,3);
    meanstrengths=mean(meanmat);
    load(nbs_effects);
    %see if there's a correlation between effect size and node strength
    [r_effectsize_centrality,p_r_effectsize_centrality]=corr(pointsz',meanstrengths')
    %compare Gtot in nodes with non-zero effects to zero-effect nodes
    sv=sort(pointsz);
    min_stren=0; %sv(round(0.9*length(sv)));
    top10=pointsz>min_stren;
    bottom90=top10==0;
    top10_stren=meanstrengths(top10);
    bottom90_stren=meanstrengths(bottom90);
    top10_pointsz=pointsz(top10);
    bottom90_pointsz=pointsz(bottom90);
    p_ranksum_effectsize_centralityu=ranksum(top10_stren,bottom90_stren)
    effect_group=zeros(1,length(pointsz));
    effect_group(bottom90)=1;
    effect_group(top10)=2;
    %plot
    plot(pointsz,meanstrengths,'ko')
    yL = get(gca,'YLim');
    line([min_stren min_stren],yL,'Color','r');
    figure();
    boxplot(meanstrengths,effect_group);
    set(gca,'FontSize',20);
    xlabel('Node Group');ylabel('Total node communicability');
    set(gca,'XTickLabel',{'No Effect','Effect'})
end