function  [stressp, stress]=extraplstress( obj, xi, eta)
% Extrapolate stress from conventional gaussian points to any point within
% the element with the local coordinates (xi, eta)
%stressp   : Effective stresses
%stress    : Total stresses
% The function should hold an column vector of xi and eta
% make sure xi and eta are column vectors.
if size(xi,2)>1
    xi=xi';
end
if size(eta,2)>1
    eta=eta';
end
% change the local coordinates from xi-eta to r-s
r = sqrt(3)*xi;
s = sqrt(3)*eta;
N = [(1-r).*(1-s),(1+r).*(1-s),(1+r).*(1+s),(1-r).*(1+s)]/4;
%% Important bug, do not add the if condition. Will stop updating the stress
% and pore pressure at the gauss points and nodes 011620
% if isempty(obj.GaussPntDictM(1).Stressp)
    u=obj.Un2i;      % elemental iterative total displacement vector, u_sup(n+1)_sub(i+1)
    p=obj.Pn2i;      % elemental iterative total pore pressure vector, p_sup(n+1)_sub(i+1)
    for ig=1:length(obj.GaussPntDictM)
        obj.GaussPntDictM(ig)=obj.GaussPntDictM(ig).matsu(u,p);
    end
% end
stressp_gauss=[obj.GaussPntDictM.Stressp];
stress_gauss=[obj.GaussPntDictM.Stress];
% Note the transpose is necessary so that
% each row represent the stresses [sx,sy,sxy,sz] at one point
% each column represent the same component from four gaussian points
stressp=N*stressp_gauss';   % Effective stresses
stress=N*stress_gauss';     % Total stresses
end

