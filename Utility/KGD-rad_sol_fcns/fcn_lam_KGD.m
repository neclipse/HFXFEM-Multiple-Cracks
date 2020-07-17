function lam = fcn_lam_KGD(Kh,Ch,alp)


lamK=0.5;%475;% for K and Kt vertex
lamM=0.588;%581;%for M vertex
lamMt=0.520;%for Mt vertex
 
%zeroth approximation for the efficiency
lam0=0.5;
del=(1+fcn_Delta_p(Kh,Ch,0))/2;
B0=@(x,p1,p2) beta(p1,p2).*(1-betainc(x,p1,p2));
fcn_B2=2.^(1+del).*(B0(1/2,lam0+1,1+del))+Ch.*alp.^(3/2).*beta(alp,3/2);
eta0=1-Ch.*alp.^(3/2).*beta(alp,3/2)./fcn_B2;

%lambda interpolation
pK=Kh.^(3);
peta=eta0;
lam=lamM*(1-pK).*peta+lamMt*(1-pK).*(1-peta)+lamK*pK;


 end

