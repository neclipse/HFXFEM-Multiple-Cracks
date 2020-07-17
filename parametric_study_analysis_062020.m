% choose what cases need to be compared
clear;clc;
addpath(genpath('.\Utility'));
cases=14:22;
% cases=[4,3,2,1];   % Tini
% cases=[5,3,6,7];   % Tkrg
% cases=[8,2,9,10];    % Lcr
% cases=[11,12,1,13];% Dmax
R=cell(size(cases));
for ic=1:length(cases)
    path='C:\Users\chuan25\Exciting PHD- Heavy stuff\TSL_Parametric_study_Task2_0602-Renumbered\';
    filename=strcat(path,'case_',num2str(cases(ic)),'.mat');
    R{ic}=load(filename);
end
cases=cases-13; % renumbered for publication 060920
%% Focus on the Final results
%% Compare the time history gg
X=zeros(length(cases),1);
L=zeros(length(cases),1);
W=zeros(length(cases),1);
P=zeros(length(cases),1);
M=zeros(length(cases),1);
CV=zeros(length(cases),1);
for ic=1:length(cases)
    % Generate plots from Results_t [NT,NL,NCMOD,NCMP,NCV,NLV];
    % length
    L(ic)=R{ic}.Results(end,2);
    W(ic)=R{ic}.Results(end,3);
    P(ic)=R{ic}.Results(end,4);
    X(ic)=cases(ic);
    M(ic)=R{ic}.GBINP.Biot_mod;
end
figurename='Impacts of shape';
h1=figure('Name',figurename,'NumberTitle','off');
% L
L_change=(L-L(1))/L(1)*100;
subplot(2,2,1)
bar(X,L_change,0.5,'r')
title('Final Crack length Change')
xlabel('Case Id')
ylabel('Comparing to case-1 (%)')
xticks(cases)
yticks(-2:2:10)
set(gca,'fontsize',13)

% % W
% subplot(2,2,2)
% plot(X,W,'--ks',...
%     'LineWidth',2,...
%     'MarkerSize',10)
% title('Final Crack Mouth Aperture (mm)','FontSize',12)
% xticks(cases)

%Net P
subplot(2,2,2)
P_change=(P-P(1))/P(1)*100;
bar(X,P_change,0.5)
title('Final Net Crack Mouth Pressure Change')
xlabel('Case Id')
ylabel('Comparing to case-1 (%)')
xticks(cases)
set(gca,'fontsize',13)

% %Crack Shape
% subplot(2,2,4)
% plot(X,P,'--gs',...
%     'LineWidth',2,...
%     'MarkerSize',10)
% title('Final Crack Volume (m)','FontSize',12)
% xticks(cases)
% Crack shape
subplot(2,2,[3,4])
for ic=1:length(cases)
    % Generate crackshape from final postprocess
    postprocess=R{ic}.postdict(end);
    snapshot(postprocess,'crackgeo',1,'basecoordinate',30);
end
hold off
set(gca,'fontsize',13)
title('Final Crack Shape')
xlabel('Distance from crack mouth (m)')
ylabel('Crack lip displacement (mm)')
legs=strcat(repmat('case-',length(cases),1),num2str(cases(:)));
legend(legs,'Location','west')
set(h1, 'Position',[48,48,830,580]);

%% exporting the final figure
export_fig 'Organized results\Parametric_study_052820\Impacts of TSL shapes_0620.tif' -m3.125 -transparent

