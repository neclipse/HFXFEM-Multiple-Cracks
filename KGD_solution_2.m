%% Analytical solution considering fluid lag
% Gather the inputs
agl=assembleglobalinputs();
t1=linspace(0,1,10);
t2=linspace(1,10,20);
t=unique([t1,t2]);
qline=0.011248593925759;              % Unit flux rate, which is volumetric rate/area: [m/s], negative for source
% the sectional area of the perforated wellbore: wellbore 3.75inch diameter (0.085725 m) ,
% 0.35 inch (0.00889m ) peforation hole; gives:
q=0.5*qline*0.00889;            % flow rate: [m^2/s] 1e-4
S=agl.inistress(2);         % inistress in y direction.
CMOD=2.36*(agl.muf*(1-agl.nu^2)*q^3/agl.E)^(1/6)*t.^(1/3);
CMOD=CMOD*1e3; % mm
L=0.539*(agl.E*q^3/agl.muf/(1-agl.nu^2))^(1/6)*t.^(2/3);
CMP=1.09*(agl.E^2*agl.muf/(1-agl.nu^2)^2)^(1/3)*t.^(-1/3)-S;
CMP=CMP*1e3; %MPa

%% Numerical solution
startpoint=2;
endpoint=45;
NT=[Step1.Postprocess(startpoint:endpoint).Inc];
NL=zeros(1,endpoint-startpoint+1);
NCMOD=NL;
NCMP=NL;
for i=startpoint:endpoint
    mouthind=length(Step1.Postprocess(i).EnrichItems{1}.IntPoints)/2+3;
    NL(i-startpoint+1)=Step1.Postprocess(i).EnrichItems{1}.Length;
    NCMOD(i-startpoint+1)=Step1.Postprocess(i).EnrichItems{1}.Aperture(mouthind);
    NCMP(i-startpoint+1)=Step1.Postprocess(i).EnrichItems{1}.Pfrack(mouthind);
end
NL=(NL-0.1)/2;  % 0.1 is the initial perforated length
NCMOD=(NCMOD-agl.minaperture)*1e3; % mm
NCMP=NCMP*1e3; %MPa

%% Comparison Plots
% length
figure()
plot(NT,NL,'*',t,L,'-','LineWidth',2);
xlabel('Time (s)')
ylabel('Crack length (m)')
legend('XFEM','KGD2')
% aperture
figure()
plot(NT,NCMOD,'*',t,CMOD,'-','LineWidth',2);
xlabel('Time (s)')
ylabel('Crack Mouth Opening Displacement (mm)')
legend('XFEM','KGD2')
%CMP
figure()
plot(NT,NCMP,'*',t,CMP,'-','LineWidth',2);
xlabel('Time (s)')
ylabel('Crack Mouth Pressure (MPa)')
legend('XFEM','KGD2')