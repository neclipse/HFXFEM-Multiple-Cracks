function  obj=matsu(obj,disp,p)
%MATSU Stress State Update method of GaussPnt_LE_UP
%   Elastic stress is simply the product of strain and elastic tangent
%   matrix
%   TENSILE STRESS IS DEFINED POSITIVE
%% Preparing       
% retrieve material properties from element class
% global Delastic inistressp Biot_alpha inipore;              % initial effective stress state, initial pore pressure, and alpha
%% poroelastic step
dirac=[1;1;0;1];
obj.Strainc=obj.Bmat*disp;
stresspc=obj.GBINP.Delastic*obj.Strainc;      % relative effective stress change comparing to the initial state
obj.Stressp=stresspc+obj.GBINP.inistressp;
obj.P=obj.Np*p+obj.GBINP.inipore;     % current total pore pressure=current excess pressure + initial total pore pressure
obj.U=obj.Nu*disp;                       
obj.Stress=obj.Stressp-obj.GBINP.Biot_alpha*dirac*obj.P;
end


 