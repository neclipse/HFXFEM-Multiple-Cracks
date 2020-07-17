function  epnotifier( obj )
%EPNOTIFIER method of Class NewRapItr
% To notify if current load increment has been marked by plastic behavior
% The notification is transferred from lower level to this level,i.e., from
% gauss point to element(crtstif), to the linear system(SparseSysCreater) and then to this
% newton-raphson iterator. At linear system, the EPFlag is set to zero when
% number of plastic elements exceeds 1% of total number of elements.
    fprintf('The %d increment is completed in %d iterations\n',obj.IInc, obj.IItr);
if obj.LinSysCrt.EPArea~=0
    fprintf('An area of %5.1e square meter has yielded after No.%d load increment\n',obj.LinSysCrt.EPArea,obj.IInc);
end

