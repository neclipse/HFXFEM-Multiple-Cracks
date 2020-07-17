function obj=assembleglobalinputs(icase)
% This fucntion is to speed up the the access of global variables
% This file is set to for "main-03152020_Abaqus_Ruhrsandstone.m"
% 4.2 Storage-toughness dominated regime
%% Parameters for parallel runs
tinis=[0.8	0.8	0.8	0.8	0.8	0.8	0.8	0.8	0.8	0.8	0.8	0.8]*1e-3;
tmaxs=[1.2	1.2	1.2	2.1	2.1	2.1	1.2	1.2	1.2	2.1	2.1	2.1]*1e-3;
lcrs =[0.4	0.4	0.4	0.4	0.4	0.4	0.4	0.4	0.4 0.4 0.4 0.4];
dcrs =[2 2	2	2	2	2	2	2	2	2	 2	2]*1e-4;
Es=[13	15	18	13	15	18	15.96	15.96	15.96	15.96	15.96	15.96];
nus=[0.219	0.219	0.219	0.219	0.219	0.219	0.16	0.25	0.28	0.16	0.25	0.28];
%% Parameters
obj=struct();
% -- Rock Properties(GN, GPa, m, s)
obj.E=Es(icase);
obj.nu=nus(icase);
% rhos=2e-6;                                % 1kg->1N, 1000kg=1e-6 GN
% rhow=1e-6;                                % Gg/m^3
% poro=0.2;                                  % porosity not the void ratio
% obj.Density=poro*rhow+(1-poro)*rhos;
obj.Density=0;
obj.lcr=lcrs(icase);
obj.dcr=dcrs(icase);  % critical crack displacement where cohesion vanishes (m) from Khoei
% (2011), which will be used as initial crack disp for perforated mode      
% obj.threshold=0.5e-3;      % set to a large value so that stationary crack is obtained
obj.threshold=tinis(icase);       % tini
obj.tmax=tmaxs(icase);           % tkrg
obj.Gc=0.5*(obj.threshold*obj.lcr+obj.tmax)*obj.dcr;        %GN.m
% lcoh=obj.E*obj.Gc/(1-obj.nu^2)/obj.tmax^2;                  % Estimated cohesive zone size based on Hillerborg et al. 
% % tmax is, however, usually size dependent according to Bazant, and Kazemi
% % 1990 and 1991
obj.lini=0;                 % Dimensionless parameter FOR initial crack separation, NOT USED BY BILINEAR NOR UNIFIED TRACTION LAW.
obj.perfaperture=1e-4;
obj.minaperture=5e-5;
obj.mu=50;                   % dynamic viscosity of the fluid, unit cp, mpa.s
obj.mul=0.1;                % dynamic viscosity of leakoff fluid,, unit cp, mpa.s
obj.poro=0.2;              
% obj.Kf=0.0138;                % Kf, bulk modulus of the fluid phase, GPa?
obj.Kf=2.1;
obj.sgmH=-9e-3;
obj.sgmh=-2e-3;
obj.sgmv=-1e-2;
obj.inipore=0;
obj.G= 0.5*obj.E/(1+obj.nu);
obj.lambda=2*obj.nu*obj.G/(1-2*obj.nu);
obj.K=obj.lambda+2*obj.G/3;
obj.Ks=34;               % Ks, bulk modulus of the solid phase, GPa 37.778
obj.Biot_alpha=1-obj.K/obj.Ks;
% obj.Biot_alpha=0;
obj.Biot_mod_crack=obj.Ks/(obj.Biot_alpha-obj.poro*(1-obj.Ks/obj.Kf));
obj.Cstar_crack=1/obj.Biot_mod_crack;    
%(Questionable biot_modulus in the literature, did not use that number on 011020)
obj.Biot_mod=obj.Biot_mod_crack;
obj.Cstar=1/obj.Biot_mod;
% Effective compressibility, (GPa^-1)
obj.Delastic=[obj.lambda+2*obj.G,obj.lambda,0,obj.lambda;
                obj.lambda,obj.lambda+2*obj.G,0,obj.lambda;
                0,0,obj.G,0; 
                obj.lambda,obj.lambda,0,obj.lambda+2*obj.G];
obj.muf=obj.mu*1e-12;       % dynamic viscosity of the fracturing fluid, unit (GPa.s)
obj.mulf=obj.mul*1e-12;       % dynamic viscosity of the leakoff fluid, unit (GPa.s)
obj.kmat=[1,0;0,1]*1e-17/obj.mulf;  % m^2/(GPa.s)
obj.kmat_crack=[1,0;0,1]*1e-17/obj.muf; % 
% obj.kmat=obj.k*9.87e-16/obj.muf;
% Initial total stress state, [sgmx,sgmy,tauxy,sgmz]
obj.inistress=[obj.sgmH;obj.sgmh;0;obj.sgmv];
dirac=[1;1;0;1];
obj.inistressp=obj.inistress+dirac*obj.Biot_alpha*obj.inipore;
end