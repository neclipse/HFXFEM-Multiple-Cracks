function [ N,N_xi,N_eta ] = shape( n, xi, eta )
%% SHAPE function
%%
% This function provides shape function vector _N_ and its first order derivative
% with respect to local corrdinates, (xi, eta), for different type of element with _n_ nodes. 
if ismember(n,[3,4,6,8])
    switch(n)
        case 3
            N=[1-xi-eta,xi,eta];
            N_xi=[-1,1,0];
            N_eta=[-1,0,1];
        case 4
            N=[(1-xi).*(1-eta),(1+xi).*(1-eta),(1+xi).*(1+eta),(1-xi).*(1+eta)]/4;
            N_xi=[eta-1,1-eta,1+eta,-1-eta]/4;
            N_eta=[xi-1,-1-xi,1+xi,1-xi]/4;
        case 6
            N=[(1-xi-eta)*(1-2*xi-2*eta),xi*(2*xi-1),eta*(2*eta-1),...
                4*(1-xi-eta)*xi,4*xi*eta,4*(1-xi-eta)*eta];
            N_xi=[-3+4*xi+4*eta,4*xi-1,0,4-8*xi-4*eta,4*eta,-4*eta];
            N_eta=[-3+4*xi+4*eta,0,4*eta-1,-4*xi,4*xi,4-8*eta-4*xi];
        case 8
            N=[-0.25*(1-xi).*(1-eta).*(1+xi+eta),0.5*(1-xi.^2).*(1-eta),...
               -0.25*(1+xi).*(1-eta).*(1-xi+eta),0.5*(1-eta.^2).*(1+xi),...
               -0.25*(1+xi).*(1+eta).*(1-xi-eta),0.5*(1-xi.^2).*(1+eta),...
               -0.25*(1-xi).*(1+eta).*(1+xi-eta),0.5*(1-eta.^2).*(1-xi)];
            N_xi=[0.25*(1-eta)*(2*xi+eta),-xi*(1-eta),0.25*(1-eta)*(2*xi-eta),...
                  0.5*(1-eta^2),0.25*(1+eta)*(2*xi+eta),-xi*(1+eta),...
                  0.25*(1+eta)*(2*xi-eta),-0.5*(1-eta^2)];
            N_eta=[0.25*(1-xi)*(2*eta+xi),-0.5*(1-xi^2),0.25*(1+xi)*(2*eta-xi),...
                   -eta*(1+xi),0.25*(1+xi)*(2*eta+xi),0.5*(1-xi^2),...
                   0.25*(1-xi)*(2*eta-xi),-eta*(1-xi)];               
    end
else
    error('The number of nodes is unacceptable.\n')
end
           
end

