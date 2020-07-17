function obj=assembleglobalinputs()
% This fucntion is to speed up the the access of global variables
obj=struct;
% -- Rock Properties(GN, GPa, m, s)
obj.E=29.120;
obj.nu=0.12;
% rhos=2.8e-6;                                % Gg/m^3, times acceleration gives MN
% rhow=1e-6;                                  % Gg/m^3
% poro=0.04;                                  % porosity not the void ratio
% rho=poro*rhow+(1-poro)*rhos;
% % Density=rho;
obj.Density=0;
% The Hillberborg cohesive model parameters were NOT SERIOUSLY CALIBRATED YET 02212019
obj.lcr=4e-1;
obj.dcr=0;  
% critical crack displacement where cohesion vanishes (m) from Khoei
% (2011), which will be used as initial crack disp for perforated mode      
obj.threshold=1;      % set to a large value so that stationary crack is obtained
obj.tmax=0;           % Dimensionless parameter to determine if the linear softening begins at the current opening 
% Gc=0.5*tmax*dcr;
% lcoh=E*Gc/(1-nu^2)/tmax^2;                  % Estimated cohesive zone size based on Hillerborg et al. 
% % tmax is, however, usually size dependent according to Bazant, and Kazemi
% % 1990 and 1991
obj.lini=0;                 % Dimensionless parameter FOR initial crack separation, NOT USED BY BILINEAR NOR UNIFIED TRACTION LAW.
% storage factor, generally positive as of filter cake and damaged crack surface with lower permeablilty
storage_factor=0;   
obj.skin=1/(1+storage_factor);
obj.mu=1;                   % dynamic viscosity of the fluid, unit cp
obj.poro=4e-2;              
obj.Kf=2.1;                 % Kf, bulk modulus of the fluid phase, GPa
obj.Ks=36;                  % Ks, bulk modulus of the solid phase, GPa
obj.k=[0.2,0;0,0.2];        % permeability in Darcy units, md
obj.sgmH=0;
obj.sgmh=0;
obj.sgmv=0;
obj.inipore=0;
obj.G= 0.5*obj.E/(1+obj.nu);
obj.lambda=2*obj.nu*obj.G/(1-2*obj.nu);
obj.K=obj.lambda+2*obj.G/3;
obj.Biot_alpha=1-obj.K/obj.Ks; 
obj.Cstar=obj.poro/obj.Kf+(obj.Biot_alpha-obj.poro)/obj.Ks;     
                            % Effective compressibility, (GPa^-1)
obj.Delastic=[obj.lambda+2*obj.G,obj.lambda,0,obj.lambda;
                obj.lambda,obj.lambda+2*obj.G,0,obj.lambda;
                0,0,obj.G,0; 
                obj.lambda,obj.lambda,0,obj.lambda+2*obj.G];
obj.muf=obj.mu*1e-12;       % dynamic viscosity of the fluid, unit (GPa.s) 
obj.kmat=obj.k*9.87e-16/obj.muf;
% Initial total stress state, [sgmx,sgmy,tauxy,sgmz]
obj.inistress=[obj.sgmH;obj.sgmh;0;obj.sgmv];
dirac=[1;1;0;1];
obj.inistressp=obj.inistress+dirac*obj.Biot_alpha*obj.inipore;
end