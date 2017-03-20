clear;clc;close all;
addpath(genpath('/home/david/Projects/matlab_scripts'));

%% USER DEFINED INPUTS
ANCOVAdesign=       load('NBSprep/designANCOVA164_123male.txt');
inmat=              'COMmats/123connectome_scale125.mat';
load('randomforests/input_files/123male_COM_scale125_NOnuisancereg.mat');

logtransform=       0; % 0 or 1
ranknormalize=      1; % if 1, overrules logtransform

regress_nuisance=   1;

% set nuisance covariates here
age=ANCOVAdesign(:,3);
TCV=ANCOVAdesign(:,4);
nuisance_predictors=[ones(length(age),1) age TCV];
% set effect of interest here (diagnosis, behavior, etc.)
effect=ANCOVAdesign(:,2);

%% BEGIN SCRIPT

mean_comm = mean(Typical,1);

% % get rid of top 10 comm values cuz they're huge
for i = 1:1
    removevec = mean_comm==max(mean_comm);
    Autism(:,removevec)=[];Typical(:,removevec)=[];mean_comm(removevec)=[];
end

mean_ASDcomm = mean(Autism,1);
std_ASDcomm = std(Autism,1);
std_comm = std(Typical,1);

[Y,I]=sort(mean_comm);
Y_std = std_comm(I);
Y_ASD_std = std_ASDcomm(I);
Y_ASD = mean_ASDcomm(I);

plot(Y,'b'); hold on; plot(Y-Y_ASD,'r'); hold off;
zero = line(xlim, [0 0],'Color',[.5 .5 .5],'LineStyle','--','LineWidth',1);


%% BEGIN SCRIPT
%load data
load(inmat);
[sx, sy, sz]=size(subject_array_3D);

%% for each subject, convert upper triangle of residuals matrix to array, and add to subject_array_2D
%% which is subjects down rows and features (i.e. connections) across columns
inds=find(triu(ones(sx,sx),1));
subject_array_2D=zeros(sz,length(inds));
for s=1:sz
    submat=subject_array_3D(:,:,s)+subject_array_3D(:,:,s)';
    subject_array_2D(s,:)=submat(inds);
end

%% loop through connections and compute ANCOVA t-stats
numcs = size(subject_array_2D,2);
tstats = zeros(numcs,1);pvals= zeros(numcs,1);nodes= zeros(numcs,2);
for i=1:size(subject_array_2D,2)
    outcome = subject_array_2D(:,i);
    
    if ranknormalize == 1
        rank = tiedrank(outcome);
        p = rank / ( length(rank) + 1 ); %# +1 to avoid Inf for the max point
        outcome = norminv( p, 0, 1 ) * std(outcome) + mean(outcome);
    end
    
    if regress_nuisance == 1
        b=nuisance_predictors\outcome;
        regressed_y=outcome-nuisance_predictors*b + mean(outcome);
        [~, P, ~, T_STATS] = ttest2(regressed_y(effect==1), regressed_y(effect==-1),'vartype','unequal');
        nodes(i,:)=round(10000*[T_STATS.tstat P])/10000;
    else
        [~, P, ~, T_STATS] = ttest2(oucome(effect==1), outcome(effect==-1),'vartype','unequal');
        nodes(i,:)=round(10000*[T_STATS.tstat P])/10000;
    end
    tstats(i)=nodes(i,1);
    pvals(i)=nodes(i,2);
end

hold on
plot(0.1*(tstats > 3.5),'ko')
plot(-0.1*(tstats < -3.5),'ko')
hold off