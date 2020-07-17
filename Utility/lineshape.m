function [ Nline,Nline_xi ] = lineshape( n, xi)
%% LINESHAPE function
%%
% This function provides line shape function vector _N_ and its first order derivative
% with respect to local corrdinates, _xi_, for different type of element with _n_ nodes. 
if ismember(n,[3,4,6,8])
switch(n)
    case 3
        Nline=[1-xi,0,xi,0;0,1-xi,0,xi];
        Nline_xi=[-1,1];
    case 4
        Nline=[1-xi,0,1+xi,0;0,1-xi,0,1+xi]/2;
        Nline_xi=[-1/2,1/2];
    case 6
        Nline=[2*(1-xi)*(1/2-xi),0,4*xi*(1-xi),0,2*xi*(xi-1/2),0;...
               0,2*(1-xi)*(1/2-xi),0,4*xi*(1-xi),0,2*xi*(xi-1/2)];
        Nline_xi=[4*xi-3,4-8*xi,4*xi-1];
    case 8
        Nline=[xi*(xi-1),0,2*(1-xi^2),0,xi*(xi+1),0;0,xi*(xi-1),0,2*(1-xi^2),0,xi*(xi+1)]/2;
        Nline_xi=[xi-1/2,-2*xi,xi+1/2];   
end
end        
end

