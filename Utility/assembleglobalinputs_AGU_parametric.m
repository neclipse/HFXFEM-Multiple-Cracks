function obj=assembleglobalinputs(varargin)
% This fucntion is to speed up the the access of global variables
% Parametric study for AGU conference
% global icase;
% iicase=icase-26; % start from 27
obj=struct;
% if ~isempty(varargin)
%     icase=varargin{1};
%     names={'kmat','skin'};

%     inputs=[kcases;skincases];
%     for ip=1:length(names)
%         switch names{ip}
%             case 'kmat'
%                 obj.k=[inputs(ip,icase),0;0,inputs(ip,icase)];
%             case 'skin'
%                 obj.skin=inputs(ip,icase);
%         end
%     end
% end
% kcases=[10,10,10,1,1,1,0.1,0.1,0.1,0.01,0.01,0.01];
% skincases=[1,0.5,0.1,1,0.5,0.1,1,0.5,0.1,1,0.5,0.1];
% tinis=[1	1	1	1	0.6	0.8	1.2	0.6	0.8	1.2	1	1]*1e-3;
% tmaxs=[1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.5	1.8]*1e-3;
% lcrs=[0.1	0.3	0.5	0.7	0.3	0.3	0.3	0.5	0.375	0.25	0.3	0.3];
% dcrs=[9.23077E-05	0.00008	7.05882E-05	6.31579E-05	8.69565E-05	8.33333E-05	7.69231E-05	0.00008	0.00008	0.00008	6.66667E-05	5.71429E-05];
% obj.threshold=tinis(iicase);
% obj.tmax=tmaxs(iicase);
% obj.lcr=lcrs(iicase);
% obj.dcr=dcrs(iicase);
% -- Rock Properties(GN, GPa, m, s)
obj.E=15.96;
obj.nu=0.12;
% rhos=2.8e-6;                                % Gg/m^3, times acceleration gives MN
% rhow=1e-6;                                  % Gg/m^3
% poro=0.04;                                  % porosity not the void ratio
% rho=poro*rhow+(1-poro)*rhos;
% % Density=rho;
obj.Density=0;
% The Hillberborg cohesive model parameters were NOT SERIOUSLY CALIBRATED YET 02212019
obj.lcr=0.3;
obj.dcr=0.00008;               % critical crack displacement where cohesion vanishes (m) from Khoei (2011) 
obj.threshold=0.5e-3;
obj.tmax=0.6e-3;
obj.Gc=0.5*obj.tmax*obj.dcr+0.5*obj.threshold*obj.lcr*obj.dcr;
obj.lcoh=obj.E*obj.Gc/(1-obj.nu^2)/obj.tmax^2;                  % Estimated cohesive zone size based on Hillerborg et al. 
% % tmax is, however, usually size dependent according to Bazant, and Kazemi
% % 1990 and 1991
obj.lini=0;                 % Dimensionless parameter FOR initial crack separation. SHOULD BE LESS THAN OBJ.LCR for 'bilinear' and 'unified' case 
obj.perfaperture=1e-5;
obj.minaperture=1e-7;
% storage factor, generally positive as of filter cake and damaged crack surface with lower permeablilty
storage_factor=-0.9;   % after 11/21/19, do not change storage factor again. Keep it zero
obj.skin=1/(1+storage_factor);
obj.mu=50;                   % dynamic viscosity of the fluid, unit cp
obj.poro=0.19;              
obj.Kf=3;                 % Kf, bulk modulus of the fluid phase, GPa
obj.Ks=36;                  % Ks, bulk modulus of the solid phase, GPa
obj.k=[0.01,0;0,0.01];      % permeability in Darcy units, md
obj.sgmH=-0.0009;
obj.sgmh=-0.0002;           % The compressive stress should be negative.
obj.sgmv=0;
obj.inipore=0;        % Atmosphere pressure
obj.G= 0.5*obj.E/(1+obj.nu);
obj.lambda=2*obj.nu*obj.G/(1-2*obj.nu);
obj.K=obj.lambda+2*obj.G/3;
obj.Biot_alpha=1-obj.K/obj.Ks;
% obj.Biot_alpha=0;
obj.Cstar=obj.poro/obj.Kf+(obj.Biot_alpha-obj.poro)/obj.Ks;     
                            % Effective compressibility, (GPa^-1)
obj.Delastic=[obj.lambda+2*obj.G,obj.lambda,0,obj.lambda;
                obj.lambda,obj.lambda+2*obj.G,0,obj.lambda;
                0,0,obj.G,0; 
                obj.lambda,obj.lambda,0,obj.lambda+2*obj.G];
obj.muf=obj.mu*1e-12;       % dynamic viscosity of the fluid, unit (GPa.s) 
obj.kmat=obj.k*9.87e-16/obj.muf;
obj.inistress=[obj.sgmH;obj.sgmh;0;obj.sgmv];
dirac=[1;1;0;1];
obj.inistressp=obj.inistress+dirac*obj.Biot_alpha*obj.inipore;
end