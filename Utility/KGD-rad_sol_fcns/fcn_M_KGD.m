function M = fcn_M_KGD(rho,s)

%elasticity kernel for a KGD crack
M=s./(s.^2-rho.^2);

end

