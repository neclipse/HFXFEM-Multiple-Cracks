function switching(obj, mode)
%SWITCHING Method of NewRapItr Class 
%   Switching between current value and last converged value 
%   according to the switching modes
%-- Mode "1": Store current values of some parameters as the last
%converged values. The parameters include U, P,Ut1, Ut2, Pt1
%This situation occurs when the convergence criterion is met, at the end 
%of the iterations for the current load increment.

%-- Mode "2": Reset the current values of some parameters to last converged value.
%Internal load vector to the last converged value;
%LHS to the last converged value;
%Call GaussPnt.Switching(2) for all gauss points.
%New External load vector and new RHS; 
%IItr=0 and Convergence flag=0. This scenario happens when the returning is
%moving to the next load increment or increment cutting is activated.

% Switching values of parameters of gaussian points class

% Switching values of parameters of NewRapItr class
obj.LinSysCrt.switching(mode);
if mode==1
    obj.DtO=obj.Dt;
else
    obj.ConvFlag=0;
    obj.DivFlag=0;
    obj.CutFlag=1;
    obj.IItr=0;
end
end
