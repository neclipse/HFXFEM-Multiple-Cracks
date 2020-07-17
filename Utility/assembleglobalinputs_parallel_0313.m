function obj=assembleglobalinputs(icase)
% This fucntion is to speed up the the access of global variables
% This file is set to for "main-01082020_Carrier.m"
% 4.2 Storage-toughness dominated regime
%% Parameters for parallel runs
tinis=[0.5, 1, 1.5, 2]*1e-3;
% tmaxs=[1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.5	1.8]*1e-3;
% lcrs=[0.1	0.3	0.5	0.7	0.3	0.3	0.3	0.5	0.375	0.25	0.3	0.3];
% dcrs=[9.23077E-05	0.00008	7.05882E-05	6.31579E-05	8.69565E-05	8.33333E-05	7.69231E-05	0.00008	0.00008	0.00008	6.66667E-05	5.71429E-05];
%% Parameters
obj=struct();
% -- Rock Properties(GN, GPa, m, s)
obj.E=17;
obj.nu=0.2;
% rhos=2e-6;                                % 1kg->1N, 1000kg=1e-6 GN
% rhow=1e-6;                                % Gg/m^3
% poro=0.2;                                  % porosity not the void ratio
% obj.Density=poro*rhow+(1-poro)*rhos;
obj.Density=0;
% The Hillberborg cohesive model parameters were NOT SERIOUSLY CALIBRATED YET 02212019
obj.lcr=0;
obj.dcr=1.92e-4;  
% critical crack displacement where cohesion vanishes (m) from Khoei
% (2011), which will be used as initial crack disp for perforated mode      
% obj.threshold=0.5e-3;      % set to a large value so that stationary crack is obtained
obj.threshold=tinis(icase);
obj.tmax=obj.threshold;           % Dimensionless parameter to determine if the linear softening begins at the current opening 
obj.Gc=0.5*obj.tmax*obj.dcr;        %
% lcoh=obj.E*obj.Gc/(1-obj.nu^2)/obj.tmax^2;                  % Estimated cohesive zone size based on Hillerborg et al. 
% % tmax is, however, usually size dependent according to Bazant, and Kazemi
% % 1990 and 1991
obj.lini=0;                 % Dimensionless parameter FOR initial crack separation, NOT USED BY BILINEAR NOR UNIFIED TRACTION LAW.
obj.perfaperture=1e-4;
obj.minaperture=1e-5;
obj.mu=100;                   % dynamic viscosity of the fluid, unit cp, mpa.s
% obj.mul=0.1;                % dynamic viscosity of leakoff fluid,, unit cp, mpa.s
obj.poro=0.2;              
% obj.Kf=0.0138;                % Kf, bulk modulus of the fluid phase, GPa?
obj.Kf=3.1;
% small Kf will lead to very small bulk modulus, i.e., pure solid-like
% material. However, one should note that the current fomrulation also use
% this Kf for the fracturing fluid, which is not contradictory to the
% assumption of "incompressible newtonian fluid" adopted by the analytical
% benchmark solutions. Therefore, for current stage, we will use
% elem.crtstif_enriched_1Dflow where compressibility is ignored and f
% factor for is set as 1 for ideal parallel plates.
% obj.k=[0.2,0;0,0.2];        % permeability in Darcy units, md
% The normal stress is consistent with literature, but the other two
% components are assigned to guarantee the crack propagate horizontally.
obj.sgmH=-9e-3;
obj.sgmh=0;
obj.sgmv=-1e-2;
obj.inipore=0;
obj.G= 0.5*obj.E/(1+obj.nu);
obj.lambda=2*obj.nu*obj.G/(1-2*obj.nu);
obj.K=obj.lambda+2*obj.G/3;
obj.Ks=37.778;               % Ks, bulk modulus of the solid phase, GPa 37.778
obj.Biot_alpha=1-obj.K/obj.Ks;
% obj.Biot_alpha=0;
obj.Biot_mod_crack=obj.Ks/(obj.Biot_alpha-obj.poro*(1-obj.Ks/obj.Kf));
obj.Cstar_crack=1/obj.Biot_mod_crack;    
%(Questionable biot_modulus in the literature, did not use that number on 011020)
obj.Biot_mod=0.0687;
obj.Cstar=1/obj.Biot_mod;
% Effective compressibility, (GPa^-1)
obj.Delastic=[obj.lambda+2*obj.G,obj.lambda,0,obj.lambda;
                obj.lambda,obj.lambda+2*obj.G,0,obj.lambda;
                0,0,obj.G,0; 
                obj.lambda,obj.lambda,0,obj.lambda+2*obj.G];
obj.muf=obj.mu*1e-12;       % dynamic viscosity of the fluid, unit (GPa.s)
% obj.mulf=obj.mul*1e-12;       % dynamic viscosity of the fluid, unit (GPa.s)
obj.kmat=[1,0;0,1]*1e-16/obj.muf;  % m^2/(GPa.s)
obj.kmat_crack=[1,0;0,1]*1e-16/obj.muf; % 
% obj.kmat=obj.k*9.87e-16/obj.muf;
% Initial total stress state, [sgmx,sgmy,tauxy,sgmz]
obj.inistress=[obj.sgmH;obj.sgmh;0;obj.sgmv];
dirac=[1;1;0;1];
obj.inistressp=obj.inistress+dirac*obj.Biot_alpha*obj.inipore;
end