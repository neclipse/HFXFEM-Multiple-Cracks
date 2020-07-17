function p=fcn_pres_cal(del,lam,rho,type)

w=(1-rho).^del.*(1+rho).^lam/2.^lam;

ds=rho(2)-rho(1);

[Rho,S]=meshgrid(rho,rho);

if type==1
    M = fcn_M_KGD(Rho,S+ds/2)-fcn_M_KGD(Rho,S-ds/2);
elseif type==2
    M = fcn_M_rad(Rho,S+ds/2)-fcn_M_rad(Rho,S-ds/2);   
end

p=M'*w/(2*pi);

%p(end)=2*p(end-1)-p(end-2);

end