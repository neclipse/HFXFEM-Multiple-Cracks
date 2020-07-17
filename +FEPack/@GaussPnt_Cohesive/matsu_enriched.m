function obj=matsu_enriched(obj,us,ue,ps,pe)
%stress update for cohesive gaussian points along crack line
%   us are the standard displacement 
%   ue are the enriched displacement
%   ps are the standard pore pressure 
%   pe are the enriched pore pressure
%% Preparing       
% retrieve material properties from element class
% global Delastic inistressp Biot_alpha inipore;

%% poroelstic step
dirac=[1;1;0;1];
obj.Strainc=obj.Bmat*us+obj.Bmatenr*ue;
stresspc=obj.GBINP.Delastic*obj.Strainc;      % relative effective stress change comparing to the initial state
obj.Stressp=stresspc+obj.GBINP.inistressp;
obj.P=obj.Np*ps+obj.Npenr*pe+obj.GBINP.inipore;     % current total pore pressure
obj.Stress=obj.Stressp-dirac*obj.P;                 % The biot alpha is not needed here. 013120
end

