function switching(obj,mode)
%SWITCHING Method of GaussPnt Class 
%   Switching between current value and last converged value 
%   according to the switching modes
%-- Mode "1": Storing current values of some parameters as the last
%converged values. The parameters include stress, strain, hardening
%parameter, and logical flags. This situation occurs when the convergence
%criterion is met, at the end of the iterations for the current load
%increment.
%-- Mode "2": Reset the current values of those parameters by the last
%converged values. This scenario happens when moving to the next load
%increment or increment cutting is activated.
if mode==1
    obj.StresscO=obj.Stressc;
    obj.ElaStnO=obj.ElaStn;
    obj.PlaStnO=obj.PlaStn;
    obj.EPBarO=obj.EPBar;
    obj.EPFlagO=obj.EPFlag;
elseif mode==2
    obj.Stressc=obj.StresscO;
    obj.ElaStn=obj.ElaStnO;
    obj.PlaStn=obj.PlaStnO;
    obj.EPBar=obj.EPBarO;
    obj.EPFlag=obj.EPFlagO;
end

