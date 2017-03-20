clear all;close all;
addpath(genpath('/home/david/Projects/matlab_scripts'));

%% I/O - matrix inputs and base directories
matdir='/media/data1/CMP/analyses/';
mtype_in={'NBSprep/subgroups/'}; %{'NBSprep/','NBSprep/exclude0s_','COMmats/'}; %,'COMmats/exclude0s'};
mtype_out={'social_poscor_COM_'}; %{'STREN_','exclude0s_STREN_','COM_'}; %,'exclude0s_COM_'};
matpref='all_maleASD_COM_ranknormalize_'; %'123connectome_';
parc='scale125'; %GordonParc_reslice or scale125
designfile='NBSprep/design_files_behavior/all_maleASD_designANCOVA_6factors_try2_social.txt'; %'NBSprep/designANCOVA164_123male.txt';

%% I/O - NBS parameters
thresh=2.6:.2:4.8; %F/t-stat thresholds in NBS
UI.method.ui='Run NBS';
UI.test.ui='t-test'; %F-test or t-test
UI.size.ui='Intensity';
UI.perms.ui='2000';
UI.alpha.ui='.9';
UI.contrast.ui='[0 0 0 1]'; %[0 1 0 0] for TYP>ASD (or M>F), [0 -1 0 0] for ASD>TYP (or F>M)
UI.design.ui=[matdir designfile];
UI.exchange.ui=[''];
UI.node_coor.ui=[''];
UI.node_label.ui=[''];

%% if you want only a subgraph of the network, set it here
%subgraph_name='somatomotor'; %only matters if you set a subgraph below
%subgraph_nodes={'CCa_L','PMCm_L','M1_L','S1_L','S2_L','PCs_L','CCa_R','PMCm_R','M1_R','S1_R','S2_R','PCs_R'}; %somatomotor sparse
subgraph_nodes={}; %use this to look across the whole matrix

%% plot settings
clrscheme={'k'}; %{'m','r','c','b'}; %must have at least as many elements as mtype_in
fontsize=20;
outstem=['/media/data1/CMP/analyses/NBSresults/' matpref parc '_' UI.test.ui '_thresh' num2str(min(thresh)) 'to' num2str(max(thresh))];

%% BEGIN SCRIPT %%
%% loop through matrix types
for m=1:length(mtype_in)
load([matdir mtype_in{m} matpref parc '.mat']);

%% get the subgraph if you specified one
if ~isempty(subgraph_nodes)
    %% load region labels
    %%%%%%%% THIS IS NOT RIGHT. YOU NEED TO INPUT AN ACTUAL ROI INFO FILE
    ROI_info_file=[matdir 'NBSprep/parcellation_info/names_' parc '.txt'];
    %%%%%%%% EITHER THAT OR CHANGE THIS SECTION TO TAKE JUST THE NAMES LIST
    fid=fopen(ROI_info_file);
    ROI_info=textscan(fid,'%d %d%d%d %s %s%s %s%s %s%s','emptyValue',0,'HeaderLines',1,'TreatAsEmpty','na');
    fclose(fid);
    labels=ROI_info{1,5};

    %% get edges of the subgraph
    [sx, sy, sz]=size(subject_array_3D);
    mod=zeros(sx,1);
    keep_conns3D=zeros(sx, sy, sz);
    for n=1:length(subgraph_nodes)
        mod(strcmpi(subgraph_nodes{n},labels))=1;
    end
    keep_conns2D=(mod*1)*(mod*1)';
    for s=1:sz; keep_conns3D(:,:,s)=keep_conns2D; end
    %apply edges to the input matrix
    subject_array_3D=subject_array_3D.*keep_conns3D;
    %rename outstem
    outstem=[outstem subgraph_name];
end

%% save your matrix as a tmp file for NBS input (deleted later)
save(['tmp' num2str(m) '.mat'],'subject_array_3D');
UI.matrices.ui=['tmp' num2str(m) '.mat']; %define input matrix

%% loop through t thresholds
c=0;
for t=thresh
    c=c+1;
    UI.thresh.ui=num2str(t); %define t-stat threshold
    disp(['T-stat threshold = ' num2str(t)]);
    NBSrun(UI); %run NBS
    global nbs;
    p=nbs.NBS.pval;
    comp(:,:,c,m)=full(nbs.NBS.con_mat{find(p==min(p))}); %get component with min p, store it
    if c > 1
        tstab(c-1,m)=1-sum(sum(abs(comp(:,:,c,m)-comp(:,:,c-1,m))))/sum(sum(comp(:,:,c-1,m))); %compute t-stability (proportion of sig. connections retained)
        thresh2(c-1)=thresh(c-1);
    end
    pval(c,m)=round(min(p)*1000)/1000; %round p to thousandth, store it
    numcons(c,m)=sum(sum(comp(:,:,c,m))); %count connections in component, store it
end
delete(['tmp' num2str(m) '.mat']);
end

%% plot # connections
figure();
for m=1:length(mtype_in)
    semilogy(thresh,numcons(:,m),clrscheme{m},'LineWidth',2);
    hold on;
end
hold off;
set(gca,'LineWidth',1,'Box','off','FontName','Arial','FontSize',fontsize);
axis([min(thresh),max(thresh),1,max(numcons(:))]);
savefig([outstem '_numconns.fig']);
print([outstem '_numconns.png'], '-dpng','-r300');

%% plot p-values
figure();
for m=1:length(mtype_in)
    semilogy(thresh,pval(:,m),clrscheme{m},'LineWidth',2);
    hold on;
end
hold off;
set(gca,'LineWidth',1,'Box','off','FontName','Arial','FontSize',fontsize);
set(gca,'YGrid','on')
zero = line(xlim, [.05 .05],'Color',[0.1 0.1 0.1],'LineStyle','--','LineWidth',1);
zero = line(xlim, [.01 .01],'Color',[0.1 0.1 0.1],'LineStyle','--','LineWidth',1);
axis([min(thresh),max(thresh),.001,.5]);
savefig([outstem '_pvals.fig']);
print([outstem '_pvals.png'],'-dpng','-r300');

%% plot t-stability
figure();
for m=1:length(mtype_in)
    plot(thresh2,tstab(:,m),clrscheme{m},'LineWidth',2);
    hold on;
end
hold off;
set(gca,'LineWidth',1,'Box','off','FontName','Arial','FontSize',fontsize);
axis([-inf,inf,0,1]);
savefig([outstem '_tstab.fig']);
print([outstem '_tstab.png'],'-dpng','-r300');

close all;