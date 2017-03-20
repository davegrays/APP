clear all;close all;
addpath(genpath('/home/david/Projects/matlab_scripts'));

%% I/O
load('randomforests/input_files/ranknormal_123male_COM_scale125.mat');
ANCOVAdesign = load('NBSprep/design_files_behavior/all_maleASD_designANCOVA_6factors_try2_social.txt');
% ANCOVAdesign = load('NBSprep/design_files_behavior/all_maleASD_designANCOVA_NVDQneg.txt');
% ANCOVAdesign = load('NBSprep/design_files_behavior/all_maleASD_designANCOVA_6factors_try2_rbeh.txt');
social = ANCOVAdesign(:,4);
mean_social = mean(social);
nSamples = length(social);endg1=60;

%% plot lambda tuning curve
cvfit = cvglmnet(Autism(1:endg1,:), social(1:endg1,:));
cvglmnetPlot(cvfit);
lambda = 0.15 %cvfit.lambda_min %could try 0.25 for social, 0.15 for rbeh, 
coefs = cvglmnetCoef(cvfit, lambda); coefs = coefs(2:end);
disp([num2str(sum(coefs~=0)) ' connections included']);
pred_social = cvglmnetPredict(cvfit, Autism(endg1+1:end,:), lambda);
disp(['correlation using ' num2str(endg1) ' for training = ' num2str(corr(pred_social, social(endg1+1:end,:)))]);

leave_how_many_out = 1;
num_perms = 3;

%% run
orig_social = social;
num_total = length(social);
num_test=leave_how_many_out;
num_train=num_total-num_test;
trainvec_start=[ones(num_train,1); zeros(num_test,1)];

for pn=1:num_perms %the very first permutation is the observed data
    disp(['permutation ' num2str(pn)]);
    if pn == 1; permvec = [1:num_total];
    else permvec=randperm(num_total); end
    social = orig_social(permvec,:);
    
    for cv=1:num_total %do loocv to get a predicted value for each subject
        trainvec_loocv = ones(num_total,1);
        trainvec_loocv(cv) = 0; trainvec_loocv = trainvec_loocv == 1;
        testvec_loocv = (1 - trainvec_loocv) == 1;
        
        %centered_Autism_train = Autism(trainvec_loocv,:) - repmat(mean(Autism(trainvec_loocv,:),1),sum(trainvec_loocv),1);
        %centered_Autism_test = Autism(testvec_loocv,:) - mean(Autism(testvec_loocv,:));
        centered_social_train = social(trainvec_loocv) - mean(social(trainvec_loocv));
        
        fit1 = glmnet(Autism(trainvec_loocv,:),centered_social_train);
        pred_social_loocv(cv,1) = glmnetPredict(fit1,Autism(testvec_loocv,:),lambda);
    end
    
    if pn == 1; r_obs = corr(pred_social_loocv,social);disp(['observed r = ' num2str(r_obs)]);
    else r_null(pn-1) = corr(pred_social_loocv,social); end
end

disp(['r = ' num2str(r_obs) ', p-value = ' num2str(mean(r_obs < r_null))]);

%% list connections involved
matdir='/media/data1/CMP/analyses/';
matpref='COMmats/ranknormal_123connectome_';
parc='scale125'; %GordonParc_reslice or scale125
load([matdir matpref parc '.mat']);
[sx, sy, sz]=size(subject_array_3D);
[x,y]=find(triu(ones(sx,sx),1));

roidir=[matdir 'NBSprep/parcellation_info/'];
ROIinfofile = [roidir 'ROI_INFO_' parc '_rescaled.txt'];
fid=fopen(ROIinfofile);
ROI_info=textscan(fid,'%d %d%d%d %s %s%s %s%s %s%s','emptyValue',0,'HeaderLines',1,'TreatAsEmpty','na');
fclose(fid);
labels=ROI_info{1,5};

disp('connections included are:')
cinds=find(coefs);
for c=1:length(cinds)
    disp([labels(x(cinds(c))) ' - ' labels(y(cinds(c)))]);
end
