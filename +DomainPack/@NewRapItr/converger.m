function  converger( obj,varargin )
%CONVERGER method of Returner class
%   Check the convergence of Newton-Raphson iteration for current load
%   increment, by calculating the residual load vector and comparing to the
%   assigned tolerance. If not converged, the residual load vector is
%   assigned to the RHS for the linear systerm creater to carry on to the
%   next iteration.
if ~isempty(varargin)
   psddofs=varargin{1}; 
else
   psddofs=[];
end
residual=obj.ResLoadVec;
external=obj.ExtLoadVec;
residual(psddofs)=0;
external(psddofs)=0;
resmax=max(abs(residual));
resnorm=norm(residual);
extnorm=norm(external);
ratio=resnorm/extnorm;
% Check convergence criteria
if ratio<obj.Tolerance||resmax<obj.Tolerance/1000
    obj.ConvFlag=1;
else
    obj.LinSysCrt.RHS=obj.ResLoadVec;
end
if obj.IItr>1
   if ratio>20*obj.Ratio || resmax>20*obj.ResMax 
      obj.DivFlag=1; 
   end
end
obj.Ratio=ratio;
obj.ResMax=resmax;
end

