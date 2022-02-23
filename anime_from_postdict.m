% Prepare the workspace
clear;clc;
caseid= 1;
% load the results
filename1=strcat('C:\Users\chuan\Google Drive\Exciting Research\Writings\Efficient HM-XFEM model with complex fracture network\Results\Enhanced rerun-0321\case_',num2str(caseid),'_newangles5.mat');
load(filename1);
%% prepare the figure canvas
% h1=figure('NumberTitle','off','Position',[40,40,640*1.5,500*1.5]);
h1=figure('NumberTitle','off','Position',[50,60,500,600]);
fs=16;
crackdict=[];
% Initialize the Frame storage place
framenum=length(postdict);
Frames(framenum)=struct('cdata',[],'colormap',[]);
%% start capturing movie frames
for iframe=framenum
    % clear current axis
    cla;
    hold on
    for ienr=1:length(postdict(1).EnrichItems)
        crackgeo=postdict(iframe).EnrichItems{ienr}.Mygeo;
        crackdict=[crackdict,crackgeo];
        if iframe==framenum
            if ienr==1 || ienr>3
                crackgeo.plotme('tipflag',true);
            else
                crackgeo.plotme;
            end
        else
            crackgeo.plotme;
        end
    end
    titlestr='Evolution of fracture network';
    title(titlestr,'FontSize',fs);
    axis('equal')
    xlabel('X axis (m)')
    ylabel('Y axis (m)')
    xlim([0,8])
    ylim([24,38])
    % Set x and y font sizes.
    ax=gca;
    ax.XAxis.FontSize = fs-2;
    ax.YAxis.FontSize = fs-2;
    legend('HF1','HF2','NF1','NF2','NF3','FontSize',fs-2,'Location','East')
    hold off;
    Frames(iframe)=getframe(h1);
end
%% save Frame as movie
vname=strcat('C:\Users\chuan\Google Drive\Exciting Research\Writings\Efficient HM-XFEM model with complex fracture network\Results\Enhanced rerun-0321\evolution_case_',num2str(caseid),'_newangles5.mp4');
v= VideoWriter(vname,'MPEG-4');
v.FrameRate=2;
open(v);
writeVideo(v,Frames)
close(v);
