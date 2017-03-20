clear all;close all;
addpath(genpath('/home/david/Projects/matlab_scripts'));

%% I/O - matrix inputs and base directories

matdir='/media/data1/CMP/analyses/';
parc='scale125'; %GordonParc_reslice or scale125
roidir=[matdir 'NBSprep/parcellation_info/'];
outdir=[matdir 'NBSresults/caret_files/'];
load('Func_Adults12_CostStrucAdults14_250_anat_mod_vec.mat');
funcmodules=Ci_ordered_subcort;
funcmodulenames=modulenames_ordered_subcort;

% mtype_in={'NBSprep/subgroups/'};%{'COMmats/'} %{'NBSprep/subgroups/'}; %{'NBSprep/','NBSprep/exclude0s_','COMmats/','COMmats/exclude0s'};
% mtype_out={'CBCi_'};%{'social_poscor_COM_loggedFX'};%{'ranknormal_COM_'} %{'rbeh_negcor_COM_'}; %{'STREN_','exclude0s_STREN_','COM_','exclude0s_COM_'};
% matpref='CBCexists_maleASD_COM_ranknormal'; %'all_maleASD_COM_ranknormalize_';%'ranknormal_123connectome_' %'all_maleASD_COM_ranknormalize_'; %'123connectome_'; 

mtype_in={'NBSprep/time2age/'};%{'NBSprep/time2age/'} %{'COMmats/'} %{'NBSprep/subgroups/'}; %{'NBSprep/','NBSprep/exclude0s_','COMmats/','COMmats/exclude0s'};
mtype_out={'tmp'};%{'social_poscor_COM_loggedFX'};%{'ranknormal_COM_'} %{'rbeh_negcor_COM_'}; %{'STREN_','exclude0s_STREN_','COM_','exclude0s_COM_'};
matpref='74male_COMscale125_ranknormal'; %'74male_COMscale125_ranknormal' %'ranknormal_123connectome_' %'all_maleASD_COM_ranknormalize_' %'all_maleASD_COM_ranknormalize_'; %'123connectome_'; 

%% set the nodes of the subgraph you want to analyze here:
%the whole matrix
subgraph_nodes={};
%scale125 cortical
%subgraph_nodes={'rh.lateralorbitofrontal_1', 'rh.lateralorbitofrontal_2', 'rh.lateralorbitofrontal_3', 'rh.lateralorbitofrontal_4', 'rh.parsorbitalis_1', 'rh.frontalpole_1', 'rh.medialorbitofrontal_1', 'rh.medialorbitofrontal_2', 'rh.medialorbitofrontal_3', 'rh.parstriangularis_1', 'rh.parstriangularis_2', 'rh.parsopercularis_1', 'rh.parsopercularis_2', 'rh.rostralmiddlefrontal_1', 'rh.rostralmiddlefrontal_2', 'rh.rostralmiddlefrontal_3', 'rh.rostralmiddlefrontal_4', 'rh.rostralmiddlefrontal_5', 'rh.rostralmiddlefrontal_6', 'rh.superiorfrontal_1', 'rh.superiorfrontal_2', 'rh.superiorfrontal_3', 'rh.superiorfrontal_4', 'rh.superiorfrontal_5', 'rh.superiorfrontal_6', 'rh.superiorfrontal_7', 'rh.superiorfrontal_8', 'rh.caudalmiddlefrontal_1', 'rh.caudalmiddlefrontal_2', 'rh.caudalmiddlefrontal_3', 'rh.precentral_1', 'rh.precentral_2', 'rh.precentral_3', 'rh.precentral_4', 'rh.precentral_5', 'rh.precentral_6', 'rh.paracentral_1', 'rh.paracentral_2', 'rh.paracentral_3', 'rh.rostralanteriorcingulate_1', 'rh.caudalanteriorcingulate_1', 'rh.posteriorcingulate_1', 'rh.posteriorcingulate_2', 'rh.isthmuscingulate_1', 'rh.postcentral_1', 'rh.postcentral_2', 'rh.postcentral_3', 'rh.postcentral_4', 'rh.postcentral_5', 'rh.supramarginal_1', 'rh.supramarginal_2', 'rh.supramarginal_3', 'rh.supramarginal_4', 'rh.superiorparietal_1', 'rh.superiorparietal_2', 'rh.superiorparietal_3', 'rh.superiorparietal_4', 'rh.superiorparietal_5', 'rh.superiorparietal_6', 'rh.superiorparietal_7', 'rh.inferiorparietal_1', 'rh.inferiorparietal_2', 'rh.inferiorparietal_3', 'rh.inferiorparietal_4', 'rh.inferiorparietal_5', 'rh.inferiorparietal_6', 'rh.precuneus_1', 'rh.precuneus_2', 'rh.precuneus_3', 'rh.precuneus_4', 'rh.precuneus_5', 'rh.cuneus_1', 'rh.cuneus_2', 'rh.pericalcarine_1', 'rh.pericalcarine_2', 'rh.lateraloccipital_1', 'rh.lateraloccipital_2', 'rh.lateraloccipital_3', 'rh.lateraloccipital_4', 'rh.lateraloccipital_5', 'rh.lingual_1', 'rh.lingual_2', 'rh.lingual_3', 'rh.fusiform_1', 'rh.fusiform_2', 'rh.fusiform_3', 'rh.fusiform_4', 'rh.parahippocampal_1', 'rh.entorhinal_1', 'rh.temporalpole_1', 'rh.inferiortemporal_1', 'rh.inferiortemporal_2', 'rh.inferiortemporal_3', 'rh.inferiortemporal_4', 'rh.middletemporal_1', 'rh.middletemporal_2', 'rh.middletemporal_3', 'rh.middletemporal_4', 'rh.bankssts_1', 'rh.superiortemporal_1', 'rh.superiortemporal_2', 'rh.superiortemporal_3', 'rh.superiortemporal_4', 'rh.superiortemporal_5', 'rh.transversetemporal_1', 'rh.insula_1', 'rh.insula_2', 'rh.insula_3', 'lh.lateralorbitofrontal_1', 'lh.lateralorbitofrontal_2', 'lh.lateralorbitofrontal_3', 'lh.lateralorbitofrontal_4', 'lh.parsorbitalis_1', 'lh.frontalpole_1', 'lh.medialorbitofrontal_1', 'lh.medialorbitofrontal_2', 'lh.parstriangularis_1', 'lh.parsopercularis_1', 'lh.parsopercularis_2', 'lh.rostralmiddlefrontal_1', 'lh.rostralmiddlefrontal_2', 'lh.rostralmiddlefrontal_3', 'lh.rostralmiddlefrontal_4', 'lh.rostralmiddlefrontal_5', 'lh.rostralmiddlefrontal_6', 'lh.superiorfrontal_1', 'lh.superiorfrontal_2', 'lh.superiorfrontal_3', 'lh.superiorfrontal_4', 'lh.superiorfrontal_5', 'lh.superiorfrontal_6', 'lh.superiorfrontal_7', 'lh.superiorfrontal_8', 'lh.superiorfrontal_9', 'lh.caudalmiddlefrontal_1', 'lh.caudalmiddlefrontal_2', 'lh.caudalmiddlefrontal_3', 'lh.precentral_1', 'lh.precentral_2', 'lh.precentral_3', 'lh.precentral_4', 'lh.precentral_5', 'lh.precentral_6', 'lh.precentral_7', 'lh.precentral_8', 'lh.paracentral_1', 'lh.paracentral_2', 'lh.rostralanteriorcingulate_1', 'lh.caudalanteriorcingulate_1', 'lh.posteriorcingulate_1', 'lh.posteriorcingulate_2', 'lh.isthmuscingulate_1', 'lh.postcentral_1', 'lh.postcentral_2', 'lh.postcentral_3', 'lh.postcentral_4', 'lh.postcentral_5', 'lh.postcentral_6', 'lh.postcentral_7', 'lh.supramarginal_1', 'lh.supramarginal_2', 'lh.supramarginal_3', 'lh.supramarginal_4', 'lh.supramarginal_5', 'lh.superiorparietal_1', 'lh.superiorparietal_2', 'lh.superiorparietal_3', 'lh.superiorparietal_4', 'lh.superiorparietal_5', 'lh.superiorparietal_6', 'lh.superiorparietal_7', 'lh.inferiorparietal_1', 'lh.inferiorparietal_2', 'lh.inferiorparietal_3', 'lh.inferiorparietal_4', 'lh.inferiorparietal_5', 'lh.precuneus_1', 'lh.precuneus_2', 'lh.precuneus_3', 'lh.precuneus_4', 'lh.precuneus_5', 'lh.cuneus_1', 'lh.pericalcarine_1', 'lh.lateraloccipital_1', 'lh.lateraloccipital_2', 'lh.lateraloccipital_3', 'lh.lateraloccipital_4', 'lh.lateraloccipital_5', 'lh.lingual_1', 'lh.lingual_2', 'lh.lingual_3', 'lh.lingual_4', 'lh.fusiform_1', 'lh.fusiform_2', 'lh.fusiform_3', 'lh.fusiform_4', 'lh.parahippocampal_1', 'lh.entorhinal_1', 'lh.temporalpole_1', 'lh.inferiortemporal_1', 'lh.inferiortemporal_2', 'lh.inferiortemporal_3', 'lh.inferiortemporal_4', 'lh.middletemporal_1', 'lh.middletemporal_2', 'lh.middletemporal_3', 'lh.middletemporal_4', 'lh.bankssts_1', 'lh.bankssts_2', 'lh.superiortemporal_1', 'lh.superiortemporal_2', 'lh.superiortemporal_3', 'lh.superiortemporal_4', 'lh.superiortemporal_5', 'lh.transversetemporal_1', 'lh.insula_1', 'lh.insula_2', 'lh.insula_3', 'lh.insula_4'};
%subcortical
%subgraph_nodes={'Right-Thalamus-Proper','Right-Caudate','Right-Putamen','Right-Pallidum','Right-Accumbens-area','Right-Hippocampus','Right-Amygdala','Left-Thalamus-Proper','Left-Caudate','Left-Putamen','Left-Pallidum','Left-Accumbens-area','Left-Hippocampus','Left-Amygdala'};
%custom subgraph
%subgraph_nodes={'Right-Amygdala','Left-Amygdala'};

subgraph_type='nodes_to_all'; %'nodes_to_all' or 'between_nodes_only'

%% NBS - parameters
thresh=[3]; %F/t-stat thresholds in NBS
UI.method.ui='Run NBS';
UI.test.ui='t-test'; %F-test or t-test
UI.size.ui='Intensity';
UI.perms.ui='2';
UI.alpha.ui='1';
UI.contrast.ui='[0 -1 0 0]'; %[0 1 0 0] for TYP>ASD (or M>F), [0 -1 0 0] for ASD>TYP (or F>M); [0 0 0 1] for positive behavior cor; [0 0 0 -1] for negative behavior cor;
UI.design.ui=[matdir 'NBSprep/time2age/designANCOVA92_74male.txt']; %'NBSprep/time2age/designANCOVA92_74male.txt' %'NBSprep/designANCOVA164_123male.txt' %'NBSprep/design_files_behavior/all_maleASD_designANCOVA_5factors_try2_ADI.txt' substitute '6factors_try2_social' for '5factors_try2_ADI' 
UI.exchange.ui=[''];
UI.node_coor.ui=[roidir 'coords_' parc '.txt'];
UI.node_label.ui=[roidir 'names_' parc '.txt'];

%% Caret visualization - parameters
make_caret_files = 0;

ROIinfofile = [roidir 'ROI_INFO_' parc '_rescaled.txt'];
ROIinfofile_int = [roidir 'ROI_INFO_int_' parc '_rescaled.txt'];
NVnodes=0; %set to >0 if you want unconnected nodes to still appear
matfiles=1; %output .mat files?
colorscheme='hold'; %'hold' if you want the first 14 colors fixed

%% BEGIN SCRIPT
%% loop through matrix types
for m=1:length(mtype_in)
%load([matdir mtype_in{m} matpref parc '.mat']);
load([matdir mtype_in{m} matpref '.mat']);
%% load region labels
fid=fopen(ROIinfofile);
ROI_info=textscan(fid,'%d %d%d%d %s %s%s %s%s %s%s','emptyValue',0,'HeaderLines',1,'TreatAsEmpty','na');
fclose(fid);
labels=ROI_info{1,5};

%% get the subgraph if you specified one
if ~isempty(subgraph_nodes)
    [sx, sy, sz]=size(subject_array_3D);
    mod=zeros(sx,1);
    keep_conns3D=zeros(sx, sy, sz);
    for n=1:length(subgraph_nodes)
        mod(strcmpi(subgraph_nodes{n},labels))=1;
    end
    %% get edges of the subgraph
    if strcmpi(subgraph_type,'nodes_to_all')
        keep_conns2D=(mod*1)*ones(1,sx);
        keep_conns2D=keep_conns2D+keep_conns2D';
    else
        keep_conns2D=(mod*1)*(mod*1)';
    end
    keep_conns2D(1:sx+1:end)=0;
    %% apply edges to the input matrix
    for s=1:sz; keep_conns3D(:,:,s)=keep_conns2D; end
    subject_array_3D=subject_array_3D.*keep_conns3D;
end
%% save your matrix as a tmp file for NBS input (deleted later)
save(['tmp' num2str(m) '.mat'],'subject_array_3D');
UI.matrices.ui=['tmp' num2str(m) '.mat']; %define input matrix

%% loop through t thresholds
for t=thresh
    UI.thresh.ui=num2str(t); %define t-stat threshold
    stem = [outdir parc '_' mtype_out{m} UI.test.ui '_thresh' UI.thresh.ui]; %prefix for output filenames
    
    %% get the stat mat
    NBSrun(UI);
    global nbs;
    statmat=nbs.NBS.test_stat;
    if length(statmat) == 234 && length(funcmodules) == 219
        statmat(227:234,:,:)=[];statmat(:,227:234,:)=[];
        statmat(109:115,:,:)=[];statmat(:,109:115,:)=[];
    end

    %% chi-square
    statmat = community_chisquared_test(statmat,funcmodules,1,t,0.5);
    funcmodulenames
    
    %% make caret foci, focicolor and vector files; use sqrt of F-values
    caret_modules=ones(length(statmat),1)*4; %go for gold! (node colors)
    if strcmp(UI.test.ui,'F-test'); statmat=sqrt(statmat);end
    if make_caret_files == 1
        caret_networkvisualizer(ROIinfofile,ROIinfofile_int,statmat,caret_modules,stem,colorscheme,NVnodes,matfiles);
    end
    
    %% make foci files using #sig connections rather than sum of tstats
    % stem=[stem '_WeightedByNumcon'];
    % binmat=conmat;
    % binmat(binmat>0)=1;
    % binaryNodes=sum(binmat);
    % caretfoci_filemaker(ROIinfofile_int,binaryNodes,modules,stem,'hold');
end

delete(['tmp' num2str(m) '.mat']);
end