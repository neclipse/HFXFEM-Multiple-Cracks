function B = fcn_B_KGD(Kh,Ch,alp)

p=0.0;%parameter for delta calculation

del=(1+fcn_Delta_p(Kh,Ch,p))/2;

lam=fcn_lam_KGD(Kh,Ch,alp);

B0=@(x,p1,p2) beta(p1,p2).*(1-betainc(x,p1,p2));

B=2.^(1+del).*(B0(1/2,lam+1,1+del))+Ch.*alp.^(3/2).*beta(alp,3/2);

end

