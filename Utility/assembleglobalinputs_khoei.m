function obj=assembleglobalinputs()
% This fucntion is to speed up the the access of global variables
% This file is set to for "main-01082020_Carrier.m"
% 4.2 Storage-toughness dominated regime
obj=struct;
% -- Rock Properties(GN, GPa, m, s)
obj.E=15.96;
obj.nu=0.2;
% rhos=2.8e-6;                                % Gg/m^3, times acceleration gives MN
% rhow=1e-6;                                  % Gg/m^3
% poro=0.04;                                  % porosity not the void ratio
% rho=poro*rhow+(1-poro)*rhos;
% % Density=rho;
obj.Density=0;
% The Hillberborg cohesive model parameters were NOT SERIOUSLY CALIBRATED YET 02212019
obj.lcr=0; 
% critical crack displacement where cohesion vanishes (m) from Khoei
% (2011), which will be used as initial crack disp for perforated mode      
obj.threshold=1e-3;      % set to a large value so that stationary crack is obtained
obj.tmax=1e-3;           % Dimensionless parameter to determine if the linear softening begins at the current opening 
obj.Gc=1.43e-7;        %
obj.dcr=2*obj.Gc/obj.tmax; 
% lcoh=obj.E*obj.Gc/(1-obj.nu^2)/obj.tmax^2;                  % Estimated cohesive zone size based on Hillerborg et al. 
% % tmax is, however, usually size dependent according to Bazant, and Kazemi
% % 1990 and 1991
obj.lini=0;                 % Dimensionless parameter FOR initial crack separation, NOT USED BY BILINEAR NOR UNIFIED TRACTION LAW.
obj.perfaperture=1e-4;
obj.minaperture=1e-4;
% storage factor, generally positive as of filter cake and damaged crack surface with lower permeablilty
storage_factor=0;   
obj.skin=1/(1+storage_factor);
obj.mu=1;                   % dynamic viscosity of the fluid, unit cp, mpa.s
obj.poro=0.19;              
% obj.Kf=0.138;                % Kf, bulk modulus of the fluid phase, GPa?
obj.Kf=3;
% small Kf will lead to very small bulk modulus, i.e., pure solid-like
% material. However, one should note that the current fomrulation also use
% this Kf for the fracturing fluid, which is not contradictory to the
% assumption of "incompressible newtonian fluid" adopted by the analytical
% benchmark solutions. Therefore, for current stage, we will use
% elem.crtstif_enriched_1Dflow where compressibility is ignored and f
% factor for is set as 1 for ideal parallel plates.
obj.Ks=36;               % Ks, bulk modulus of the solid phase, GPa
% obj.k=[0.2,0;0,0.2];        % permeability in Darcy units, md
% The normal stress is consistent with literature, but the other two
% components are assigned to guarantee the crack propagate horizontally.
obj.sgmH=0;
obj.sgmh=0;
obj.sgmv=0;
obj.inipore=0;
obj.G= 0.5*obj.E/(1+obj.nu);
obj.lambda=2*obj.nu*obj.G/(1-2*obj.nu);
obj.K=obj.lambda+2*obj.G/3;
% obj.Biot_alpha=1-obj.K/obj.Ks;
obj.Biot_alpha=0.79;
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
obj.muf=obj.mu*1e-12;       % dynamic viscosity of the fluid, unit (GPa.s)
obj.kmat=[6,0;0,6]*1e-15/obj.muf;
obj.kmat_crack=[6,0;0,6]*1e-15/obj.muf; % only used by blending elements
% obj.kmat=obj.k*9.87e-16/obj.muf;
% Initial total stress state, [sgmx,sgmy,tauxy,sgmz]
obj.inistress=[obj.sgmH;obj.sgmh;0;obj.sgmv];
dirac=[1;1;0;1];
obj.inistressp=obj.inistress+dirac*obj.Biot_alpha*obj.inipore;
end