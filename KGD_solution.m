%% Analytical solution
% Gather the inputs
agl=assembleglobalinputs();
t=linspace(0,14,30);
% t2=linspace(1,10,20);
% t=unique([t1,t2]);
qline=0.011248593925759;              % Unit flux rate, which is volumetric rate/area: [m/s], negative for source
% the sectional area of the perforated wellbore: wellbore 3.75inch diameter (0.085725 m) ,
% 0.35 inch (0.00889m ) peforation hole; gives:
q=0.0005;            % flow rate: [m^2/s] 1e-4
S=agl.inistress(2);         % inistress in y direction.
CMOD=1.32*(agl.muf*(1-agl.nu)*q^3/agl.G)^(1/6)*t.^(1/3);
CMOD=CMOD*1e3; % mm
L=0.48*(agl.G*q^3/agl.muf/(1-agl.nu))^(1/6)*t.^(2/3);
CMP=0.955*(agl.G^3*q*agl.muf/(1-agl.nu)^3./L.^2).^(1/4)-S;
CMP=CMP*1e3; %MPa

%% Numerical solution
startpoint=1;
endpoint=length(Step1.Postprocess);
% endpoint=10;
NT=[Step1.Postprocess(startpoint:endpoint).Inc];
NL=zeros(1,endpoint-startpoint+1);
NCMOD=NL;
NCMP=NL;
LKF=NL;
CKF=NL;
for i=startpoint:endpoint
    mouthind=1;
    NL(i-startpoint+1)=Step1.Postprocess(i).EnrichItems{1}.Length;
    NCMOD(i-startpoint+1)=Step1.Postprocess(i).EnrichItems{1}.Aperture(mouthind);
    NCMP(i-startpoint+1)=Step1.Postprocess(i).EnrichItems{1}.Pfrack(mouthind);
    LKF(i-startpoint+1)=Step1.Postprocess(i).EnrichItems{1}.LeakoffVolume;
    CKF(i-startpoint+1)=Step1.Postprocess(i).EnrichItems{1}.CrackVolume;
end
NL=(NL-0.05);  % 0.05 is the initial perforated length
NCMOD=(NCMOD-agl.perfaperture)*1e3; % mm
NCMP=NCMP*1e3-3.7; %MPa

%% Comparison Plots
% length
figure()
plot(NT,NL,'*',t,L,'-','LineWidth',2);
xlabel('Time (s)')
ylabel('Crack length (m)')
legend('XFEM','KGD','Location','Northwest')
% aperture
figure()
plot(NT,NCMOD,'*',t,CMOD,'-','LineWidth',2);
xlabel('Time (s)')
ylabel('Crack Mouth Opening Displacement(mm)')
legend('XFEM','KGD','Location','Northwest')
%CMP
figure()
plot(NT,NCMP,'*',t,CMP,'-','LineWidth',2);
xlabel('Time (s)')
ylabel('Crack Mouth Pressure (MPa)')
legend('XFEM','KGD','Location','Northwest')
% leakoff 
figure()
plot(NT,LKF,'*',NT,CKF,'o')
xlabel('Sqrt of Time (s^{-1/2})')
ylabel('Leak-off rate (m^2/s)')