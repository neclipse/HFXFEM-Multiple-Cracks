function M = fcn_M(rho,s)

M=0*s;

ind1=find(rho>s);
ind2=find(rho<s);

s1=s(ind1);
rho1=rho(ind1);
s2=s(ind2);
rho2=rho(ind2);

[K1,E1]=ellipke(s1.^2./rho1.^2);
[~,E2]=ellipke(rho2.^2./s2.^2);

M(ind1)=1./rho1.*K1+rho1./(s1.^2-rho1.^2).*E1;

M(ind2)=s2./(s2.^2-rho2.^2).*E2;


end

