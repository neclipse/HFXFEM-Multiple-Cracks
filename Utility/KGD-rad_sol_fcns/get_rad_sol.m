function [R,w,p,rho,eta]=get_rad_sol(Ep,mup,Kp,Cp,Q0,t,N,plotfig)

%     Kp(Kp<1e4)=1e4;
%     mup(mup<1e-10)=1e-10;
%     Ep(Ep<1e5)=1e5;
%     Cp(Cp<1e-20)=1e-20;
%     Q0(Q0<1e-10)=1e-10;
%     t(t<1e-10)=1e-10;
    
    tmk=(mup^5*Ep^(13)*Q0^3/Kp^(18))^(1/2);

    %dimensionless parameters
    tau=t/tmk;
    phi=mup^3*Ep^(11)*Cp^4*Q0/Kp^(14);

    %determine length and width
    tau2=[tau/4 tau/2 tau];
    [Om,gamma,eff,del,lam] = rad_HF_appr(tau2,phi);

    %scales
    Lst=(Q0^3*Ep*tmk^4/mup)^(1/9);
    Eps=(mup/(Ep*tmk))^(1/3);

    rho0=linspace(0,1,N+1);
    rho=(rho0(1:end-1)+rho0(2:end))/2;
    
    %unscaled results
    R=gamma(3)*Lst;
    w=Om(3)*Eps*Lst*(1-rho).^(del(3)).*(1+rho).^(lam(3));
    p=Eps*Ep*2.^lam(3)*Om(3)./gamma(3)*fcn_pres_cal(del(3),lam(3),rho',2);  
    eta=eff(3);
    
    %determine regime of propagation
    %disp(tau);
    %disp(phi);
    %plotfig=1;
    
    if plotfig==1
        %axis limits
        taumin=-10;
        taumax=20;
        phimin=-30;
        phimax=20;
        
        %boundaries
        tmk0=4.54*1e-2;
        tmk1=2.59*1e6;

        tmmt0=7.41*1e-6;
        tmmt1=7.20*1e2;

        tkkt0=5.96*1e-8;
        tkkt1=4.81*1e2;

        tmtkt0=4.18;
        tmtkt1=2.01*1e11;
        
        %M zone
        %Taum=min(tmk0*ones(size(phi)),tmmt0./phi.^(9/14));

        %K zone
        %Phik=(tkkt0./tau).^(6/5);
        %iKK=find(tau>=tmk1);

        %Mt zone
        %Phimt=max((tmmt1./tau).^(14/9),(tau/tmtkt0).^2);

        %Kt zone
        %Taukt=max(tmtkt1*phi.^(1/2),tkkt1./phi.^(5/6));
        
        figure;
        hold on;
        
 
        taumt=(tmmt1^(14/9)*tmtkt0^2)^(9/32);
        phimt=(taumt/tmtkt0)^2;
        taukt=(tkkt1^(6/5)*tmtkt1^2)^(5/16);
        phikt=(taukt/tmtkt1)^2;
        
        %plot edge limits
        plot3(log10([tmk0 tmk0 tmmt0/(10^phimax)^(9/14)]),log10([10^phimin (tmk0/tmmt0)^(-14/9) 10^phimax]),1e5*ones(1,3),'b-');%M
        plot3(log10([tmmt1/(10^phimax)^(9/14) taumt tmtkt0*10^(phimax/2)]),log10([10^phimax phimt 10^phimax]),1e5*ones(1,3),'g-');%Mt
        plot3(log10([tmk1  tmk1 tkkt0/(10.^phimin)^(5/6)]),log10([10^phimin (tkkt0/tmk1)^(6/5) 10^phimin]),1e5*ones(3),'r-');%K
        plot3(log10([10^taumax taukt 10^taumax]),log10([(tkkt1/10^(taumax))^(6/5) phikt (10^taumax/tmtkt1)^2]),1e5*ones(1,3),'m-');%Kt

        %add labels
        text(-7.5,-16,'$M$','fontsize',30);
        text(9,-27,'$K$','fontsize',30);
        text(0,13,'$\tilde M$','fontsize',30);
        text(14,-5,'$\tilde K$','fontsize',30);  
        
        title(['$\tau=$' num2str(tau,3) ' $\phi=$' num2str(phi,3)]);

        %show the location on the phase diagram
        tau(tau<10.^taumin)=10.^taumin;
        tau(tau>10.^taumax)=10.^taumax;
        phi(phi<10.^phimin)=10.^phimin;
        phi(phi>10.^phimax)=10.^phimax;
     
        plot3(log10(tau),log10(phi),1e5*ones(1,3),'ko','markerfacecolor','k');
        
        xlabel('$\log_{10}(\tau)$');
        ylabel('$\log_{10}(\phi)$');

        
        axis([taumin taumax phimin phimax]);
        view(0,90);
    end
end

