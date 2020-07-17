function  matct( obj )
%MATCT Method of GaussPnt_DP to calculate the tangent operator consistent
%with the last converged stress returing
%   This method will use elastic trial strain and dgama for last returning
%   mapping as inputs, and output the Tangent property of the object
%% Initializing and retriving material properties
global  G K xi eta etabar H Delastic; 
Is=diag([1 1 0.5 1]);                           % The symmetric identity tensor in array form
I=[1 1 0 1];                                    % The identity tensor in array form
etds=zeros(4,1);                                % elastic trial deviatoric strain tensor in array form
dgama=obj.DGama;                                
etrial=obj.ETrial;                              % elastic trial strain tensor in array form
epflag=obj.EPFlag;
apexflag=obj.ApexFlag;
%% Tangent consistent with the apex
if epflag
    if apexflag
        alpha=xi/etabar;
        beta=xi/eta;
        afact=K*(1-K/(K+alpha*beta*H));
        Dmat=afact*(I'*I);
%% Tangent consistent with the smooth cone
    else
% -set up elastic trial deviatoric strain (physical)
        etv=etrial(1)+etrial(2)+etrial(4);
        etds([1,2,4])=etrial([1,2,4])-etv/3;
        etds(3)=etrial(3)/2;
        etdsnorm=sqrt(etds(1)^2+etds(2)^2+etds(4)^2+2*etds(3)^2);
% -unit deviatoric flow vector
        if etdsnorm
            etdsuni=etds/etdsnorm;              % unit vector of etds
        else
            etdsuni=zeros(4,1);      
        end
% -Assemble tangent operator
        aconst=1/(G+K*eta*etabar+xi^2*H);   % 'A' constant
        afact=2*G*(1-dgama/(sqrt(2)*etdsnorm));
        bfact=2*G*(dgama/(sqrt(2)*etdsnorm)-G*aconst);
        cfact=-sqrt(2)*G*K*aconst;
        dfact=K*(1-K*eta*etabar*aconst);
    % deviatoric tensor in array form
    Id=Is-(I'*I)/3;
    Dmat=afact*Id+bfact*(etdsuni*etdsuni')+cfact*(eta*(etdsuni*I)+...
        etabar*(I'*etdsuni'))+dfact*(I'*I); % WAY LARGER THAN THE D ELASTIC
    end
else
%% Elasticity matrix
Dmat=Delastic;
end
obj.Tangent=Dmat;
end

