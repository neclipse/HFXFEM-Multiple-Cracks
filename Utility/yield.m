function  yield( stress )
%YIELD Check the yield condition 
global  xi eta epbar c0 H
ps=(stress(1)+stress(2)+stress(4))/3;   % Initial pressure stress
ds([1,2,4])=stress([1,2,4])-ps;           % Initial deviatoric stress
ds(3)=stress(3);                         
% Compute elastic trial stress J2 invariant
SQRJ2=sqrt(ds(3)^2+0.5*(ds(1)^2+ds(2)^2+ds(4)^2));
cohe=c0+H*epbar;
% Check plastic admissity
phi=SQRJ2+eta*ps-xi*cohe;                     % Drucker-Prager yield surface function
if  cohe~=0
    residual=phi/abs(cohe);
end
if residual>1e-10
    error('The initial stress state is already out of pure elastic region.')
end
end

