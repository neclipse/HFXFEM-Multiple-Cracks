function B = fcn_B_rad(Kh,Ch,alp)

p=0.0;%parameter for delta calculation

del=(1+fcn_Delta_p(Kh,Ch,p))/2;

lam=fcn_lam_rad(Kh,Ch,alp);

B0=@(x,p1,p2) beta(p1,p2).*(1-betainc(x,p1,p2));

B=2.^(1+del).*(-B0(1/2,lam+1,2+del)+B0(1/2,lam+2,1+del))+Ch.*alp.^(3/2).*beta(2*alp,3/2);

end

