addpath(genpath('/remote_home/grayson/scripts/matlab_scripts/'));
clear all;close all;

%% I/O
roidir='/media/data1/CMP/analyses/NBSprep/parcellation_info/';
infolder='NBSprep';
outfolder='COMmats';
% how many subs in the 3D conmat?
snum=164;

% what parcellation do you want? 'GordonParc_reslice', 'GordonParc_NOisolates', ...
% 'scale33', 'scale60', 'scale125', etc.
parc='scale33';

% set the nodes of the subgraph you want to analyze here:
%scale125 cortical
%subgraph_nodes={'rh.lateralorbitofrontal_1', 'rh.lateralorbitofrontal_2', 'rh.lateralorbitofrontal_3', 'rh.lateralorbitofrontal_4', 'rh.parsorbitalis_1', 'rh.frontalpole_1', 'rh.medialorbitofrontal_1', 'rh.medialorbitofrontal_2', 'rh.medialorbitofrontal_3', 'rh.parstriangularis_1', 'rh.parstriangularis_2', 'rh.parsopercularis_1', 'rh.parsopercularis_2', 'rh.rostralmiddlefrontal_1', 'rh.rostralmiddlefrontal_2', 'rh.rostralmiddlefrontal_3', 'rh.rostralmiddlefrontal_4', 'rh.rostralmiddlefrontal_5', 'rh.rostralmiddlefrontal_6', 'rh.superiorfrontal_1', 'rh.superiorfrontal_2', 'rh.superiorfrontal_3', 'rh.superiorfrontal_4', 'rh.superiorfrontal_5', 'rh.superiorfrontal_6', 'rh.superiorfrontal_7', 'rh.superiorfrontal_8', 'rh.caudalmiddlefrontal_1', 'rh.caudalmiddlefrontal_2', 'rh.caudalmiddlefrontal_3', 'rh.precentral_1', 'rh.precentral_2', 'rh.precentral_3', 'rh.precentral_4', 'rh.precentral_5', 'rh.precentral_6', 'rh.paracentral_1', 'rh.paracentral_2', 'rh.paracentral_3', 'rh.rostralanteriorcingulate_1', 'rh.caudalanteriorcingulate_1', 'rh.posteriorcingulate_1', 'rh.posteriorcingulate_2', 'rh.isthmuscingulate_1', 'rh.postcentral_1', 'rh.postcentral_2', 'rh.postcentral_3', 'rh.postcentral_4', 'rh.postcentral_5', 'rh.supramarginal_1', 'rh.supramarginal_2', 'rh.supramarginal_3', 'rh.supramarginal_4', 'rh.superiorparietal_1', 'rh.superiorparietal_2', 'rh.superiorparietal_3', 'rh.superiorparietal_4', 'rh.superiorparietal_5', 'rh.superiorparietal_6', 'rh.superiorparietal_7', 'rh.inferiorparietal_1', 'rh.inferiorparietal_2', 'rh.inferiorparietal_3', 'rh.inferiorparietal_4', 'rh.inferiorparietal_5', 'rh.inferiorparietal_6', 'rh.precuneus_1', 'rh.precuneus_2', 'rh.precuneus_3', 'rh.precuneus_4', 'rh.precuneus_5', 'rh.cuneus_1', 'rh.cuneus_2', 'rh.pericalcarine_1', 'rh.pericalcarine_2', 'rh.lateraloccipital_1', 'rh.lateraloccipital_2', 'rh.lateraloccipital_3', 'rh.lateraloccipital_4', 'rh.lateraloccipital_5', 'rh.lingual_1', 'rh.lingual_2', 'rh.lingual_3', 'rh.fusiform_1', 'rh.fusiform_2', 'rh.fusiform_3', 'rh.fusiform_4', 'rh.parahippocampal_1', 'rh.entorhinal_1', 'rh.temporalpole_1', 'rh.inferiortemporal_1', 'rh.inferiortemporal_2', 'rh.inferiortemporal_3', 'rh.inferiortemporal_4', 'rh.middletemporal_1', 'rh.middletemporal_2', 'rh.middletemporal_3', 'rh.middletemporal_4', 'rh.bankssts_1', 'rh.superiortemporal_1', 'rh.superiortemporal_2', 'rh.superiortemporal_3', 'rh.superiortemporal_4', 'rh.superiortemporal_5', 'rh.transversetemporal_1', 'rh.insula_1', 'rh.insula_2', 'rh.insula_3', 'lh.lateralorbitofrontal_1', 'lh.lateralorbitofrontal_2', 'lh.lateralorbitofrontal_3', 'lh.lateralorbitofrontal_4', 'lh.parsorbitalis_1', 'lh.frontalpole_1', 'lh.medialorbitofrontal_1', 'lh.medialorbitofrontal_2', 'lh.parstriangularis_1', 'lh.parsopercularis_1', 'lh.parsopercularis_2', 'lh.rostralmiddlefrontal_1', 'lh.rostralmiddlefrontal_2', 'lh.rostralmiddlefrontal_3', 'lh.rostralmiddlefrontal_4', 'lh.rostralmiddlefrontal_5', 'lh.rostralmiddlefrontal_6', 'lh.superiorfrontal_1', 'lh.superiorfrontal_2', 'lh.superiorfrontal_3', 'lh.superiorfrontal_4', 'lh.superiorfrontal_5', 'lh.superiorfrontal_6', 'lh.superiorfrontal_7', 'lh.superiorfrontal_8', 'lh.superiorfrontal_9', 'lh.caudalmiddlefrontal_1', 'lh.caudalmiddlefrontal_2', 'lh.caudalmiddlefrontal_3', 'lh.precentral_1', 'lh.precentral_2', 'lh.precentral_3', 'lh.precentral_4', 'lh.precentral_5', 'lh.precentral_6', 'lh.precentral_7', 'lh.precentral_8', 'lh.paracentral_1', 'lh.paracentral_2', 'lh.rostralanteriorcingulate_1', 'lh.caudalanteriorcingulate_1', 'lh.posteriorcingulate_1', 'lh.posteriorcingulate_2', 'lh.isthmuscingulate_1', 'lh.postcentral_1', 'lh.postcentral_2', 'lh.postcentral_3', 'lh.postcentral_4', 'lh.postcentral_5', 'lh.postcentral_6', 'lh.postcentral_7', 'lh.supramarginal_1', 'lh.supramarginal_2', 'lh.supramarginal_3', 'lh.supramarginal_4', 'lh.supramarginal_5', 'lh.superiorparietal_1', 'lh.superiorparietal_2', 'lh.superiorparietal_3', 'lh.superiorparietal_4', 'lh.superiorparietal_5', 'lh.superiorparietal_6', 'lh.superiorparietal_7', 'lh.inferiorparietal_1', 'lh.inferiorparietal_2', 'lh.inferiorparietal_3', 'lh.inferiorparietal_4', 'lh.inferiorparietal_5', 'lh.precuneus_1', 'lh.precuneus_2', 'lh.precuneus_3', 'lh.precuneus_4', 'lh.precuneus_5', 'lh.cuneus_1', 'lh.pericalcarine_1', 'lh.lateraloccipital_1', 'lh.lateraloccipital_2', 'lh.lateraloccipital_3', 'lh.lateraloccipital_4', 'lh.lateraloccipital_5', 'lh.lingual_1', 'lh.lingual_2', 'lh.lingual_3', 'lh.lingual_4', 'lh.fusiform_1', 'lh.fusiform_2', 'lh.fusiform_3', 'lh.fusiform_4', 'lh.parahippocampal_1', 'lh.entorhinal_1', 'lh.temporalpole_1', 'lh.inferiortemporal_1', 'lh.inferiortemporal_2', 'lh.inferiortemporal_3', 'lh.inferiortemporal_4', 'lh.middletemporal_1', 'lh.middletemporal_2', 'lh.middletemporal_3', 'lh.middletemporal_4', 'lh.bankssts_1', 'lh.bankssts_2', 'lh.superiortemporal_1', 'lh.superiortemporal_2', 'lh.superiortemporal_3', 'lh.superiortemporal_4', 'lh.superiortemporal_5', 'lh.transversetemporal_1', 'lh.insula_1', 'lh.insula_2', 'lh.insula_3', 'lh.insula_4'};
%scale60 cortical
%subgraph_nodes={'rh.lateralorbitofrontal_1', 'rh.lateralorbitofrontal_2', 'rh.parsorbitalis_1', 'rh.frontalpole_1', 'rh.medialorbitofrontal_1', 'rh.medialorbitofrontal_2', 'rh.parstriangularis_1', 'rh.parsopercularis_1', 'rh.rostralmiddlefrontal_1', 'rh.rostralmiddlefrontal_2', 'rh.superiorfrontal_1', 'rh.superiorfrontal_2', 'rh.superiorfrontal_3', 'rh.superiorfrontal_4', 'rh.caudalmiddlefrontal_1', 'rh.precentral_1', 'rh.precentral_2', 'rh.precentral_3', 'rh.paracentral_1', 'rh.rostralanteriorcingulate_1', 'rh.caudalanteriorcingulate_1', 'rh.posteriorcingulate_1', 'rh.isthmuscingulate_1', 'rh.postcentral_1', 'rh.postcentral_2', 'rh.supramarginal_1', 'rh.supramarginal_2', 'rh.superiorparietal_1', 'rh.superiorparietal_2', 'rh.superiorparietal_3', 'rh.inferiorparietal_1', 'rh.inferiorparietal_2', 'rh.inferiorparietal_3', 'rh.precuneus_1', 'rh.precuneus_2', 'rh.cuneus_1', 'rh.pericalcarine_1', 'rh.lateraloccipital_1', 'rh.lateraloccipital_2', 'rh.lateraloccipital_3', 'rh.lingual_1', 'rh.lingual_2', 'rh.fusiform_1', 'rh.fusiform_2', 'rh.parahippocampal_1', 'rh.entorhinal_1', 'rh.temporalpole_1', 'rh.inferiortemporal_1', 'rh.inferiortemporal_2', 'rh.middletemporal_1', 'rh.middletemporal_2', 'rh.bankssts_1', 'rh.superiortemporal_1', 'rh.superiortemporal_2', 'rh.transversetemporal_1', 'rh.insula_1', 'rh.insula_2', 'lh.lateralorbitofrontal_1', 'lh.lateralorbitofrontal_2', 'lh.parsorbitalis_1', 'lh.frontalpole_1', 'lh.medialorbitofrontal_1', 'lh.parstriangularis_1', 'lh.parsopercularis_1', 'lh.rostralmiddlefrontal_1', 'lh.rostralmiddlefrontal_2', 'lh.rostralmiddlefrontal_3', 'lh.superiorfrontal_1', 'lh.superiorfrontal_2', 'lh.superiorfrontal_3', 'lh.superiorfrontal_4', 'lh.caudalmiddlefrontal_1', 'lh.precentral_1', 'lh.precentral_2', 'lh.precentral_3', 'lh.precentral_4', 'lh.paracentral_1', 'lh.rostralanteriorcingulate_1', 'lh.caudalanteriorcingulate_1', 'lh.posteriorcingulate_1', 'lh.isthmuscingulate_1', 'lh.postcentral_1', 'lh.postcentral_2', 'lh.postcentral_3', 'lh.supramarginal_1', 'lh.supramarginal_2', 'lh.superiorparietal_1', 'lh.superiorparietal_2', 'lh.superiorparietal_3', 'lh.inferiorparietal_1', 'lh.inferiorparietal_2', 'lh.precuneus_1', 'lh.precuneus_2', 'lh.cuneus_1', 'lh.pericalcarine_1', 'lh.lateraloccipital_1', 'lh.lateraloccipital_2', 'lh.lingual_1', 'lh.lingual_2', 'lh.fusiform_1', 'lh.fusiform_2', 'lh.parahippocampal_1', 'lh.entorhinal_1', 'lh.temporalpole_1', 'lh.inferiortemporal_1', 'lh.inferiortemporal_2', 'lh.middletemporal_1', 'lh.middletemporal_2', 'lh.bankssts_1', 'lh.superiortemporal_1', 'lh.superiortemporal_2', 'lh.transversetemporal_1', 'lh.insula_1', 'lh.insula_2'};
%scale33 cortical
%subgraph_nodes={'rh.lateralorbitofrontal', 'rh.parsorbitalis', 'rh.frontalpole', 'rh.medialorbitofrontal', 'rh.parstriangularis', 'rh.parsopercularis', 'rh.rostralmiddlefrontal', 'rh.superiorfrontal', 'rh.caudalmiddlefrontal', 'rh.precentral', 'rh.paracentral', 'rh.rostralanteriorcingulate', 'rh.caudalanteriorcingulate', 'rh.posteriorcingulate', 'rh.isthmuscingulate', 'rh.postcentral', 'rh.supramarginal', 'rh.superiorparietal', 'rh.inferiorparietal', 'rh.precuneus', 'rh.cuneus', 'rh.pericalcarine', 'rh.lateraloccipital', 'rh.lingual', 'rh.fusiform', 'rh.parahippocampal', 'rh.entorhinal', 'rh.temporalpole', 'rh.inferiortemporal', 'rh.middletemporal', 'rh.bankssts', 'rh.superiortemporal', 'rh.transversetemporal', 'rh.insula', 'lh.lateralorbitofrontal', 'lh.parsorbitalis', 'lh.frontalpole', 'lh.medialorbitofrontal', 'lh.parstriangularis', 'lh.parsopercularis', 'lh.rostralmiddlefrontal', 'lh.superiorfrontal', 'lh.caudalmiddlefrontal', 'lh.precentral', 'lh.paracentral', 'lh.rostralanteriorcingulate', 'lh.caudalanteriorcingulate', 'lh.posteriorcingulate', 'lh.isthmuscingulate', 'lh.postcentral', 'lh.supramarginal', 'lh.superiorparietal', 'lh.inferiorparietal', 'lh.precuneus', 'lh.cuneus', 'lh.pericalcarine', 'lh.lateraloccipital', 'lh.lingual', 'lh.fusiform', 'lh.parahippocampal', 'lh.entorhinal', 'lh.temporalpole', 'lh.inferiortemporal', 'lh.middletemporal', 'lh.bankssts', 'lh.superiortemporal', 'lh.transversetemporal', 'lh.insula'};
subgraph_nodes={}; %use this to look across the whole matrix

%% GO!!
cfile=['connectome_' parc '.mat'];

%% load the 3d connection mat
con3d=[infolder '/' num2str(snum) cfile];
load(con3d);

%% get the subgraph if you specified one
if ~isempty(subgraph_nodes)
    cfile=['cortical_' cfile];
    %% load region labels
    ROIinfofile = [roidir 'ROI_INFO_' parc '.txt'];
    fid=fopen(ROIinfofile);
    ROI_info=textscan(fid,'%d %d%d%d %s %s%s %s%s %s%s','emptyValue',0,'HeaderLines',1,'TreatAsEmpty','na');
    fclose(fid);
    labels=ROI_info{1,5};
    %% get node indices of the subgraph
    [sx, sy, sz]=size(subject_array_3D);
    mod=zeros(sx,1);
    for n=1:length(subgraph_nodes)
        mod(strcmpi(subgraph_nodes{n},labels))=1;mod=mod==1;
    end
    %% create new matrix (truncated in first 2 dims)
    for s=1:sz %loop through subjects
        D = subject_array_3D(:,:,s);
        g1=0;
        g1_within_mat = zeros(sum(mod),sum(mod));
        for i=1:sx %loop through nodes of original 2d matrix and build the subgraph
            if mod(i) == 1
                g1=g1+1;
                g1_within_mat(g1,:)=D(i,mod); %add a row of the subgraph
            end
        end
        new_array_3D(:,:,s) = g1_within_mat;
    end
    %reassign subject_array_3D
    subject_array_3D=new_array_3D;
end

%% convert to communicability and save
subject_array_3D=convertMAT3D_com_wei(subject_array_3D);
mkdir(outfolder)
savefile=[outfolder '/' num2str(snum) cfile];
save(savefile,'subject_array_3D');

%% ranknormalize and save
subject_array_3D = transform_data(subject_array_3D,0,1);
savefile=[outfolder '/ranknormal_' num2str(snum) cfile];
save(savefile,'subject_array_3D');

%% load the 3d connection mat, excluding 0's
con3d=[infolder '/exclude0s_' num2str(snum) cfile];
load(con3d);

%% convert to communicability and save
subject_array_3D=convertMAT3D_com_wei(subject_array_3D);
savefile=[outfolder '/exclude0s' num2str(snum) cfile];
save(savefile,'subject_array_3D');