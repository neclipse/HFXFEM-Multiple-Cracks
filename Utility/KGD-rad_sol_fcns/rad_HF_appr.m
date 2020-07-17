function [Om,gamma,eff,del,lam] = rad_HF_appr(tau,phi)
%Note that lh is the scaled radius
phi(phi<1e-30)=1e-30;
th=tau*(2^6*phi^(3/2));
Qh=8/pi*phi;

%initial guess (for viscous regime)
alp=4/9;
Kh0=0*th+1/2;
Ch0=0*th+1/2;

%alpha iteration 
for ialp=1:3
    
Res=1;
ittmax=100;
itt=0;
tol=1e-5;

while (itt<ittmax)&&(Res>tol)

    if (itt==ittmax-1)&&(ialp==3)
       disp('No convergence, wl_radHF_appr'); 
       disp(Res);
    end
    itt=itt+1;
    Kh=Kh0;
    Ch=Ch0;
    
    %dKh=1e-5;
    %dCh=1e-10;
    
    ittK=0;
    ResK=1;
    while (ittK<ittmax)&&(ResK>tol)
        ittK=ittK+1;
        fg=fcn_g_del(Kh,Ch);
        %fgK=(fcn_g_del(Kh+dKh,Ch)-fg)/dKh;
        fgK=(-fg)./(1-Kh+eps);%secant
        
        f1=Kh.^6-alp.^(1/2)./th.^(1/2).*Ch.^(3).*fg;
        f1K=6*Kh.^5-alp.^(1/2)./th.^(1/2).*Ch.^(3).*fgK;

        Kh=0.0*Kh+1.0*(Kh-f1./f1K);
        Kh(Kh<0)=1e-5;
        Kh(Kh>1)=1-1e-5;
        ResK=max(abs(f1./f1K));
    end
    

    ittC=0;
    ResC=1;
    while (ittC<ittmax)&&(ResC>tol)
       ittC=ittC+1;
       Chtest=Ch;
       
       %if max(Ch)>1
       %  Ch=th.^(3/8).*Kh./alp.^(5/8)./Qh.^(1/4).*(fcn_B(Kh,Ch,alp)./Ch).^(1/4);                     
       %else
         Ch=th.^(3/10).*Kh.^(4/5)./alp.^(1/2)./Qh.^(1/5).*(fcn_B_rad(Kh,Ch,alp)).^(1/5);          
       %end
       
       ResC=max(abs(Ch-Chtest));
    end

    Res=max(((Kh-Kh0).^2+(Ch-Ch0).^2).^(1/2));
    

    Kh0=Kh;
    Ch0=Ch;
end

sh=fcn_g_del(Kh,Ch);

%calculate length
lh=Ch.^4.*sh.^2./Kh.^(10);

%update alpha
alp=0*lh;
alp(2:end)=(log(lh(2:end))-log(lh(1:end-1)))./(log(th(2:end))-log(th(1:end-1)));
alp(1)=alp(2);


end

p=0.0;%parameter for delta calculations
del=(1+fcn_Delta_p(Kh,Ch,p))/2;

%efficiency
eff=1-Ch.*alp.^(3/2).*beta(2*alp,3/2)./fcn_B_rad(Kh,Ch,alp);

%width at the wellbore
lam=fcn_lam_rad(Kh,Ch,alp);
wha=Ch.^2.*sh./Kh.^6./2.^(lam);

%converting to original scaling
gamma=lh/(2^4*phi);
Om=wha/(2^2*phi^(1/2));


end


