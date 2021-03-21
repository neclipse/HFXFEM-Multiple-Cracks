% retrieve mesh for all
% mesh=postdict(1).EnrichItems{1}.Mygeo.Mesh;
% mesh.plotmesh;
clear;clc
%% specify the compared cases
caseids=[7,8];
figurename='Impacts of Ambient Fracture';
h1=figure('Name',figurename,'NumberTitle','off','Position',[40,40,640*1.5,500*1.5]);
%% first case
filename1=strcat('C:\Users\chuan\Google Drive\Exciting Research\Writings\Efficient HM-XFEM model with complex fracture network\Results\Enhanced run-0318\case_',num2str(caseids(1)),'.mat');
load(filename1);
subplot(1,length(caseids),1)
hold on;
crackdict=[];
% specify the ind
ind=length(postdict);
for ienr=1:length(postdict(1).EnrichItems)
    crackgeo=postdict(ind).EnrichItems{ienr}.Mygeo;
    crackdict=[crackdict,crackgeo];
    if ienr==1 || ienr==4
        crackgeo.plotme('tipflag',true);
    else
        crackgeo.plotme; 
    end
end
ax=gca;
fs=16;
% titlestr='\Lambda_{cr} = 0.2';
% titlestr='\sigma_{h} = 2MPa';
titlestr='with NF1';
% titlestr=['With NF1 at angle = 80' char(176)];
title(titlestr,'FontSize',fs);
axis('equal')
xlabel('X axis (m)')
ylabel('Y axis (m)')
xlim([0,8])
ylim([24,38])
% Set x and y font sizes.
ax.XAxis.FontSize = fs-2;
ax.YAxis.FontSize = fs-2;
legend('HF1','HF2','NF1','NF2','NF3','FontSize',fs-2,'Location','East')
% retrieve the crackdict
hold off;

%% second case
filename2=strcat('C:\Users\chuan\Google Drive\Exciting Research\Writings\Efficient HM-XFEM model with complex fracture network\Results\Enhanced run-0318\case_',num2str(caseids(2)),'.mat');
load(filename2);
% mesh=postdict(1).EnrichItems{1}.Mygeo.Mesh;
% mesh.plotmesh;
subplot(1,length(caseids),2)
hold on;
crackdict=[];
% specify the ind
ind=length(postdict);
for ienr=1:length(postdict(1).EnrichItems)
    crackgeo=postdict(ind).EnrichItems{ienr}.Mygeo;
    crackdict=[crackdict,crackgeo];
    if ienr==1 || ienr==4
        crackgeo.plotme('tipflag',true);
    else
        crackgeo.plotme; 
    end
end
ax=gca;
fs=16;
% titlestr='\Lambda_{cr} = 0.35';
% titlestr='\sigma_{h} = 3MPa';
titlestr='Without NF1';
% titlestr=['With NF1 at angle = 30' char(176)];
title(titlestr,'FontSize',fs);
axis('equal')
xlabel('X axis (m)')
ylabel('Y axis (m)')
xlim([0,8])
ylim([24,38])
% Set x and y font sizes.
ax.XAxis.FontSize = fs-2;
ax.YAxis.FontSize = fs-2;
legend('HF1','HF2','NF2','NF3','FontSize',fs-2,'Location','East')
% retrieve the crackdict
hold off;
%%
% export_fig 'C:\Users\chuan\Google Drive\Exciting Research\Writings\Efficient HM-XFEM model with complex fracture network\Results\Enhanced run-0318\Impacts of Ambient Fracture.tif' -m3.125 -transparent