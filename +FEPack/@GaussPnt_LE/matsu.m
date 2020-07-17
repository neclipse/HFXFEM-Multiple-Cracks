function  matsu(obj)
%MATSU Stress State Update method of GaussPnt_LE
%   Elastic stress is simply the product of strain and elastic tangent
%   matrix
%   TENSILE STRESS IS DEFINED POSITIVE
%% Preparing       
% retrieve material properties from element class
% global   Delastic inistress;             % initial stress state

%% Elastic step
% Last returned total elastic strain tensor + strain change at this iteration
obj.ElaStn=obj.ElaStn+obj.ItrStn;              
stress=obj.GBINP.Delastic*obj.ElaStn;
stressc=stress-obj.GBINP.inistress;
obj.Stressc=stressc;
end

 