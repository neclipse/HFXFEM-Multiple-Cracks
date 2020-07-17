function [l,w,p,xi,eta]=get_KGD_sol(Ep,mup,Kp,Cp,Q0,t,N,plotfig)
    
%     Kp(Kp<1e4)=1e4;
%     mup(mup<1e-10)=1e-10;
%     Ep(Ep<1e5)=1e5;
%     Cp(Cp<1e-20)=1e-20;
%     Q0(Q0<1e-10)=1e-10;
%     t(t<1e-10)=1e-10;
    
    tmmt=mup*Q0^3/(Ep*Cp^6);

    %dimensionless parameters
    tau=t/tmmt;
    Km=(Kp^4/(mup*Q0*Ep^3))^(1/4);


    %determine length and width
    tau2=[tau/4 tau/2 tau];
    [Om,gamma,eff,del,lam,~] = KGD_HF_appr(tau2,Km);
    
    %scales
    Lst=(mup*Q0^5/(Ep*Cp^8))^(1/2);
    Eps=Cp^2/Q0;

    xi0=linspace(0,1,N+1);
    xi=(xi0(1:end-1)+xi0(2:end))/2;

    
    %unscaled results
    l=gamma(3)*Lst;
    w=Om(3)*Eps*Lst*(1-xi).^(del(3)).*(1+xi).^(lam(3));
    p=Eps*Ep*2.^lam(3)*Om(3)./gamma(3)*fcn_pres_cal(del(3),lam(3),xi',1);  
    eta=eff(3);
    
    % to only get the crack mouth values
    w=w(1);
    p=p(1);


    %determine regime of propagation
    %disp(tau);
    %disp(Km);
    
    %plotfig=1;
    
    if plotfig==1
        %axis limits
        taumin=-30;
        taumax=25;
        Kmin=-2.5;
        Kmax=3;

        %boundaries
        tmmt0=1.21e-13;
        tmmt1=2.36e6;

        Kmk0=0.70;
        Kmk1=4.80;

        tkkt0=1.25e-14;
        tkkt1=1.76e5;

        Kmtkt0=0.90;
        Kmtkt1=4.80;

        figure;
        hold on;

        %plot edge limits
        plot3(log10([tmmt0 tmmt0 10^taumin]),log10([10^Kmin Kmk0 Kmk0]),1e5*ones(1,3),'b-');%M
        plot3(log10([tmmt1 tmmt1 10^taumax]),log10([10^Kmin Kmtkt0 Kmtkt0]),1e5*ones(1,3),'g-');%Mt
        plot3(log10([10^taumin  tkkt0*Kmk1^(4) tkkt0*1e2^(4)]),log10([Kmk1 Kmk1 10^Kmax]),1e5*ones(3),'r-');%K
        plot3(log10([10^taumax tkkt1*Kmtkt1^(4) tkkt1*1e2^(4)]),log10([Kmtkt1 Kmtkt1 10^Kmax]),1e5*ones(1,3),'m-');%Kt

        %add labels
        text(-23.5,-1.35,'$M$','fontsize',30,'Interpreter','latex');
        text(-20.5,1.85,'$K$','fontsize',30,'Interpreter','latex');
        text(14,-1.35,'$\tilde M$','fontsize',30,'Interpreter','latex');
        text(15.5,1.7,'$\tilde K$','fontsize',30,'Interpreter','latex'); 
        title(['$\tau=$' num2str(tau,3) ' $K_m=$' num2str(Km,3)]);

        %show the location on the phase diagram
        tau(tau<10.^taumin)=10.^taumin;
        tau(tau>10.^taumax)=10.^taumax;
        Km(Km<10.^Kmin)=10.^Kmin;
        Km(Km>10.^Kmax)=10.^Kmax;
     
        plot3(log10(tau),log10(Km),1e5*ones(1,3),'ko','markerfacecolor','k');

        xlabel('$$\log_{10}(\tau)$$');
        ylabel('$$\log_{10}(K_m)$$');

        axis([taumin taumax Kmin Kmax]);
        view(0,90);
    end
end
