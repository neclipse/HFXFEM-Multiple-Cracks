classdef GlobalInput < handle
    % This is a constant class to store all common constants
    properties (Constant)
        E=29.120;
        nu=0.12;
        Density=0;
        lcr=4e-1;
        dcr=5e-5;
        tmax=3.56e-3;
        lini=0;
        skin=8e-1;
        threshold=3.56e-3;
        mu=1;          % dynamic viscosity of the fluid, unit cp
        poro=4e-2;
        Kf=2.1;
        Ks=36;
        k=[0.2,0;0,0.2];
        sgmH=-1e-2;
        sgmh=-6e-3;
        sgmv=-1.1e-2;
        inipore=5e-3;
    end
    properties(Dependent)
        G
        lambda
        K
        Biot_alpha
        Cstar
        Delastic
        muf; % dynamic viscosity of the fluid, unit (GPa.s) 
        kmat
        inistress 
        inistressp
    end
    
    methods

        function obj = GlobalInput()
            
        end
        function G=get.G(obj)
           G= 0.5*obj.E/(1+obj.nu);
        end
        function lambda=get.lambda(obj)
           lambda=2*obj.nu*obj.G/(1-2*obj.nu);
        end
        function K=get.K(obj)
           K=obj.lambda+2*obj.G/3;
        end
        function Biot_alpha=get.Biot_alpha(obj)
           Biot_alpha=1-obj.K/obj.Ks; 
        end
        function Cstar=get.Cstar(obj)
           Cstar=obj.poro/obj.Kf+(obj.Biot_alpha-obj.poro)/obj.Ks;
        end
        function Delastic=get.Delastic(obj)
           Delastic=[obj.lambda+2*obj.G,obj.lambda,0,obj.lambda;
                obj.lambda,obj.lambda+2*obj.G,0,obj.lambda;
                0,0,obj.G,0; 
                obj.lambda,obj.lambda,0,obj.lambda+2*obj.G];
        end
        function muf=get.muf(obj)
           muf=obj.mu*1e-12;
        end
        function kmat=get.kmat(obj)
           kmat=obj.k*9.87e-16/obj.muf;
        end
        function inistress=get.inistress(obj)
           inistress=[obj.sgmH;obj.sgmh;0;obj.sgmv];
        end
        function inistressp=get.inistressp(obj)
           dirac=[1;1;0;1];
           inistressp=obj.inistress+dirac*obj.Biot_alpha*obj.inipore;
        end
    end
en

