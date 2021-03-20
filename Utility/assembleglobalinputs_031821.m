function obj=assembleglobalinputs(icase)
% This fucntion is to speed up the the access of global variables
% This file is set to for "main-03152020_Abaqus_Ruhrsandstone.m"
% 4.2 Storage-toughness dominated regime
%% Parameters for parallel runs
tinis=[0	0.8	0.8	0.8	0.8	0.8	0.8	0.8	0.8	0.8 0.8 0.8]*1e-3;
tkrgs=[0	1.5	1.5	1.5 2.1	2.1	2.1 2.1 1.5	1.5	1.5	1.5]*1e-3;
lcrs =[0.2	0.4	0.4	0.4	0.4	0.4	0.4	0.4	0.4 0.4 0.4 0.4];
dcrs =[1 2	2 2 2 2 2 2 2.66 2.66 2.66 2.66]*1e-4;
% Es=[15.96	15.96	15.96	15.96	15.96	15.96 15.96	15.96 15.96	15.96 15.96	15.96];
% nus=[0.219	0.219	0.219	0.219	0.219	0.219	0.219	0.219 0.219	0.219 0.219	0.219];
% Kss=[22	26	30 34 22 26	30 34 22 26	30 34]; % case 47-58
% Kfs=[1.0 1.6 2.1 2.6 1.0 1.6 2.1 2.6 1.0 1.6 2.1 2.6]; % case 59-70
poros=[0.2,0.2,0.25,0.3,0.15,0.2,0.25,0.3,0.15,0.2,0.25,0.3];% case 71-82
%% Parameters
obj=struct();
% -- Rock Properties(GN, GPa, m, s)
obj.E=15.96;
obj.nu=0.219;
% rhos=2e-6;                                % 1kg->1N, 1000kg=1e-6 GN
% rhow=1e-6;                                % Gg/m^3
% poro=0.2;                                  % porosity not the void ratio
% obj.Density=poro*rhow+(1-poro)*rhos;
obj.Density=0;
obj.lcr=lcrs(icase);
obj.dcr=dcrs(icase);  % critical crack displacement where cohesion vanishes (m) from Khoei
obj.threshold=tinis(icase);       % this is actually tini.
obj.threshold_formaxps=1e-3; % only used for maxps grow check, when I intentionally put tini=tkrg=0 for brittle material.
obj.threshold_smeared=2e-4;    % threshold to convert "smeared" crack to open.
obj.tkrg=tkrgs(icase);           % tkrg
obj.Gc=0.5*(obj.threshold*obj.lcr+obj.tkrg)*obj.dcr;        %GN.m
% lcoh=obj.E*obj.Gc/(1-obj.nu^2)/obj.tkrg^2;                  % Estimated cohesive zone size based on Hillerborg et al. 
% % tmax is, however, usually size dependent according to Bazant, and Kazemi
% % 1990 and 1991
obj.lini=0;                 % Dimensionless parameter FOR initial crack separation, NOT USED BY BILINEAR NOR UNIFIED TRACTION LAW.
obj.perfaperture=2e-4;
obj.minaperture=8e-5;
obj.mu=10;                   % dynamic viscosity of the fluid, unit cp, mpa.s
obj.mul=1;                % dynamic viscosity of leakoff fluid,, unit cp, mpa.s
obj.poro=poros(icase);              
% obj.Kf=0.0138;                % Kf, bulk modulus of the fluid phase, GPa?
obj.Kf=2.1;
obj.sgmH=-5e-3;
obj.sgmh=-2e-3;
obj.sgmv=-8e-3;
obj.inipore=0;              % total initial pore pressure. 
obj.G= 0.5*obj.E/(1+obj.nu);
obj.lambda=2*obj.nu*obj.G/(1-2*obj.nu);
obj.K=obj.lambda+2*obj.G/3;
obj.Ks=34;               % Ks, bulk modulus of the solid phase, GPa 37.778
obj.Biot_alpha=1-obj.K/obj.Ks;
% obj.Biot_alpha=0;  
% note this biot_mod should be 0.0687 when comparing to the single solid phase
obj.Biot_mod=obj.Ks/(obj.Biot_alpha-obj.poro*(1-obj.Ks/obj.Kf));
% obj.Biot_mod=0.0687;
obj.Cstar=1/obj.Biot_mod;
% obj.Biot_mod_crack=obj.Biot_mod;
obj.Biot_mod_crack=obj.Biot_mod;
obj.Cstar_crack=1/obj.Biot_mod_crack;  
% Effective compressibility, (GPa^-1)
obj.Delastic=[obj.lambda+2*obj.G,obj.lambda,0,obj.lambda;
                obj.lambda,obj.lambda+2*obj.G,0,obj.lambda;
                0,0,obj.G,0; 
                obj.lambda,obj.lambda,0,obj.lambda+2*obj.G];
obj.muf=obj.mu*1e-12;       % dynamic viscosity of the fracturing fluid, unit (GPa.s)
obj.mulf=obj.mul*1e-12;       % dynamic viscosity of the leakoff fluid, unit (GPa.s)
obj.kmat=[1,0;0,1]*1e-19/obj.mulf;  % m^2/(GPa.s)
% kmat_crack is also for the element domain flow, mulf should be used. (not significant?)
obj.kmat_crack=[1,0;0,1]*1e-17/obj.mulf; % 
% obj.kmat=obj.k*9.87e-16/obj.muf;
% Initial total stress state, [sgmx,sgmy,tauxy,sgmz]
obj.inistress=[obj.sgmH;obj.sgmh;0;obj.sgmv];
dirac=[1;1;0;1];
obj.inistressp=obj.inistress+dirac*obj.Biot_alpha*obj.inipore;
end