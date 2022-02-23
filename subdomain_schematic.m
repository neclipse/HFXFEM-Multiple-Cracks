% Example element to generate subdomain with mutliple cracks
clear; clc;
% add the home folder to the matlab search path. Already saved for future
% session, so I commented this line.
% addpath('C:\Users\chuan\Documents\HFXFEM-Verified-Github\HFXFEM-Verified-Singlecrack');
% Import class packages 
import ToolPack.*;                          % Complimentary tool classes, like Preprocessor
import FEPack.*;                            % Elementary level classes, like element class, node class and gausspnt class
import DomainPack.*;                        % Domain level classes, like Domain, LinSysCrt, NewRapItr
import EnrichPack.*;
% Step 1: generate an instance of FEPack.Elem_2d_UP, parent element
nodedict.X=[-1,1,1,-1];
nodedict.Y=[-1,-1,1,1];
nodelist=[1,2,3,4];
% crt gaussian point
p=2;
nnodes=4;
[Q,W]=gausstable(p,'QUAD');
map=[1,3,4,2];
gauss=[Q,W];
gauss=gauss(map,:);
numgauss=size(gauss,1); % number of gauss points
GBINP=assembleglobalinputs(1);
gaussdictm(1,4)=FEPack.GaussPnt_LE_UP(); %generate an void object array
for igauss=1:4
    gaussdictm(igauss)=FEPack.GaussPnt_LE_UP(gauss(igauss,1),gauss(igauss,2),gauss(igauss,3),nnodes,GBINP);
end
elem=FEPack.Elem_2d_UP(1,4,4,nodelist,nodedict,gaussdictm);
% Step 2: set enrich and assign intersections
initialmode=1;
elem.setenrich(1,initialmode); % id , initialmode
elem.setenrich(2,initialmode);
section = 4; 
% 1: one small triangle and 1 pentagon, with intersection
% 2: one small triangle and 1 pentagon, without intersection
% 3: two similar quadrilaterals, with intersection; 
% 4: two similar quadrilaterals, without intersection; 
switch section
    case 1
        localint1=[-1,-0.5;0.3,1];
        localint2=[-0.4,1;0.5,-1];
    case 2
        localint1=[-1,-0.5;0.3,1];
        localint2=[-0.3,-1;1,0.4];
    case 3
        localint1=[-0.5,-1;0.3,1];
        localint2=[-1,0;1,0];
    case 4
        localint1=[-0.5,-1;-0.3,1];
        localint2=[0.2,1;0.4,-1];
end
elem.LocalInt{1}=localint1;
elem.LocalInt{2}=localint2;
elem.GlobalInt{1}=localint1;
elem.GlobalInt{2}=localint2;
elem.linegauss(1,'unified',initialmode);
elem.linegauss(2,'unified',initialmode);
% Step 3: plotting
figure('Units','inches','OuterPosition',[2,2,5.5,5.5])
hold on
% Highlight the two crack line in the element
plot(localint1(:,1),localint1(:,2),'k','LineWidth',2.5);
plot(localint2(:,1),localint2(:,2),'k','LineWidth',2.5);
xlim([-1,1])
ylim([-1,1])
axis tight
% plot the triangulations
elem.subdomain;
% plot the integration points on domain
plot([elem.EnrichGauss.X],[elem.EnrichGauss.Y],'ob', 'MarkerSize',10,'MarkerFaceColor',[0.6,0.6,1])
% plot the line gaussian points
plot([elem.LineGaussDict{1}.X],[elem.LineGaussDict{1}.Y],'pr','MarkerSize',12, 'MarkerFaceColor',[1 .2 .2])
plot([elem.LineGaussDict{2}.X],[elem.LineGaussDict{2}.Y],'pr','MarkerSize',12, 'MarkerFaceColor',[1 .2 .2])
% Remove ticks
set(gca,'XTick',[]);
set(gca,'YTick',[]);
% 