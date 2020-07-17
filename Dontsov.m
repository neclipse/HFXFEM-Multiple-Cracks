%% from current model run
% agl=assembleglobalinputs();
%% from saved Step1
% mesh=Step1.Preprocess.Mesher;
% elemdict=Step1.ElemDict;
% agl=Step1.GBINP;
% postdict=Step1.Postprocess;
%% from saved GBINP and postdict
agl=GBINP;
%% Numerical solution
startpoint=1;
endpoint=length(postdict);
% endpoint=61;
NT=transpose([postdict(startpoint:endpoint).Inc]);
NL=zeros(size(NT));
NCMOD=NL;
NCMP=NL;
CKV=NL;
LKV=NL;
% samplenodes=[445,4836,1615,1635,1645,1705,1755,1865,1915,1945]; % mesh
% Dontsov 0403
% samplenodes=[9333,9317,9297,9277,9257,9237,9217,9207]; % after 02/24/20
% samplenodes=[9337,9317,9297,9277,9257,9237,9217];   % prior to 02/24/20
% LKF=zeros(length(NL),length(samplenodes));
% CL=zeros(1,length(samplenodes));
for i=startpoint:endpoint
    mouthind=1;
    NL(i-startpoint+1)=postdict(i).EnrichItems{1}.Length;
    NCMOD(i-startpoint+1)=postdict(i).EnrichItems{1}.Aperture(mouthind);
    NCMP(i-startpoint+1)=postdict(i).EnrichItems{1}.Pfrack(mouthind);
    LKV(i-startpoint+1)=postdict(i).EnrichItems{1}.LeakoffVolume;
    CKV(i-startpoint+1)=postdict(i).EnrichItems{1}.CrackVolume;
%     LKF(i-startpoint+1,:)=postdict(i).EnrichItems{1}.LeakoffFlux;
end
NL=(NL-0.05);  % 0.05 is the initial perforated length
NCMOD=(NCMOD-agl.perfaperture)*1e3; % mm
NCMP=NCMP*1e3; %MPa
Results_t=[NT,NL,NCMOD,NCMP,CKV,LKV];
%% Do the linear regression based leak-off rate vs time
% nLKF=cell(1,size(LKF,2));
% for ip=1:size(LKF,2)
%     elems=mesh.nodetoelem(samplenodes(ip));
%     elemdict(elems(1)).callength;
%     elength=0.5*elemdict(elems(1)).Length(1);% the horizontal length of element
%     % multipled by 0.5 means the leakoff from this node only represent half crack face
%     temp=LKF(:,ip)/elength;      % transform from surface rate m^2/s to line rate m/s
%     ind0=find(temp,1);
%     if ~isempty(ind0)
%         t0=NT(ind0);
%         dt=NT-t0;
%         x=1./sqrt(dt(ind0+5:end));  % remove the first several points as the leakoff in the begining can be extremely high
%         y=temp(ind0+5:end);
%         nLKF{ip}=[x,y];
%         CL(ip)=x\y;% 2CL actually
%     end
% end
% CL=CL(CL>=0);
%% Approximate analytical solution
% %solution for KGD and radial HF in the dimensional form
% %do not set zero values to any of the parameters, put a a very small value instead
% %material parameters
% 
% Ep=agl.E/(1-agl.nu^2)*1e9;%Pa
% mup=12*agl.mu*1e-3;%Pa*s
% Kp=4*sqrt(2/pi*agl.Gc*1e9*Ep);%Pa*m^(1/2)
% % Cp=sum(CL)/length(CL); % 2CL actually
% Cp=2*1e-6;%m/s^(1/2)
% 
% %evaluation time
% tsmall=transpose(linspace(0.01,2,15));
% tbig=transpose(linspace(2.05,30,35));%s
% t=[tsmall;tbig];
% % t=200;
% %total injection rate, not half
% Q0=0.001;%m^3/s
% H=1;%m Q=Q0/H for KGD
% 
% %number of points
% N=100;
% 
% %fracture geometry
% type=1;%1 - KGD, 2 - radial
% plotfig=0;%1 - plot parametric space figure, 0 - do not plot
% 
% 
% L=zeros(size(t));
% CMOD=L;
% CMP=L;
% %global solution
% for i=1:length(t)
%     % approximate solution
%     [L(i),CMOD(i),CMP(i),xi,eta]=get_KGD_sol(Ep,mup,Kp,Cp,Q0/H,t(i),N,plotfig);
%     % vertex solutions
% %     ind=3 ;%1 - M, 2 - Mt, 3 - K, 4 - Kt
% %     [lv,wv,pv,xiv,etav]=KGD_vert_sol(Ep,mup,Kp,Cp,Q0/H,t(i),N);
% %     L(i)=lv(ind); CMOD(i)=wv(ind); CMP(i)=pv(ind);
% end
% CMP=CMP/1e6;
% CMOD=CMOD*1000;
%% Comparison Plots
% % length
% figure()
% plot(NT,NL,'*',t,L,'-','LineWidth',2);
% xlabel('Time (s)')
% ylabel('Crack length (m)')
% legend('XFEM','Analytical','Location','Northwest')
% % aperture
% figure()
% plot(NT,NCMOD,'*',t,CMOD,'-','LineWidth',2);
% xlabel('Time (s)')
% ylabel('Crack Mouth Opening Displacement(mm)')
% legend('XFEM','Analytical','Location','Northwest')
% %CMP
% figure()
% plot(NT,NCMP,'*',t,CMP,'-','LineWidth',2);
% xlabel('Time (s)')
% ylabel('Net Crack Mouth Pressure (MPa)')
% legend('XFEM','Analytical','Location','Northwest')
% % leakoff volume against crack volume
% figure()
% plot(NT,LKV,'*',NT,CKV,'o')
% xlabel('Time (s))')
% ylabel('Leak-off volume vs. Crack Volume')
% legend('Leak-off volume','Crack volume','Location','Northwest')
% % % regression of carter's leakoff coefficient 
% figure()
% scatter(nLKF{1}(:,1),nLKF{1}(:,2))
% hold on
% plot(nLKF{1}(:,1),CL(1)*nLKF{1}(:,1))
% xlabel('{(t-t_{0})}^{-1/2}')
% for ip=2:size(LKF,2)
%     plot(1./sqrt(NT),LKF(:,ip))
%     plot
% end