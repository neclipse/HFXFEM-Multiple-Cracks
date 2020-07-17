%solution for KGD and radial HF in the dimensional form
% clear

%add folder to path
% addpath('KGD-rad_sol_fcns');

%do not set zero values to any of the parameters, put a a very small value instead
%material parameters
agl=assembleglobalinputs();
Ep=agl.E/(1-agl.nu^2)*1e9;%Pa
mup=12*agl.mu*1e-3;%Pa*s
Kp=4*sqrt(2/pi*agl.Gc*1e9*Ep);%Pa*m^(1/2)
Cp=2*1.47e-5;%m/s^(1/2)

%evaluation time
t=14;%s

%injection rate
Q0=0.001;%m^3/s
H=1;%m Q=Q0/H for KGD

%number of points
N=1000;

%fracture geometry
type=1;%1 - KGD, 2 - radial
plotfig=1;%1 - plot parametric space figure, 0 - do not plot

if type==1
    
    %global solution
    [l,w,p,xi,eta]=get_KGD_sol(Ep,mup,Kp,Cp,Q0/H,t,N,plotfig);
    
    %vertex solutions
    [lv,wv,pv,xiv,etav]=KGD_vert_sol(Ep,mup,Kp,Cp,Q0/H,t,N);
    
    ind=3;%1 - M, 2 - Mt, 3 - K, 4 - Kt
    if ind==1
        col='b';
    elseif ind==2
        col='g';
    elseif ind==3
        col='r';
    elseif ind==4
        col='m';
    end
    
    figure;
    plot(l*xi,w,'k-');
    hold on;
    plot(lv(ind)*xiv,wv(:,ind),'--','color',col);
    
    xlabel('$x$ [m]');
    ylabel('$w$ [m]');

    figure;
    hold on;
    plot(l*xi,p,'k-');
    plot(lv(ind)*xiv,pv(:,ind),'--','color',col);
    
    xlabel('$x$ [m]');
    ylabel('$p$ [Pa]');
end

if type==2
       
    %global solution
    [R,w,p,rho,eta]=get_rad_sol(Ep,mup,Kp,Cp,Q0,t,N,plotfig);
        
    %vertex solutions
    [Rv,wv,pv,rhov,etav]=rad_vert_sol(Ep,mup,Kp,Cp,Q0,t,N);
    
    ind=1;%1 - M, 2 - Mt, 3 - K, 4 - Kt
    if ind==1
        col='b';
    elseif ind==2
        col='g';
    elseif ind==3
        col='r';
    elseif ind==4
        col='m';
    end
    
    figure;
    plot(R*rho,w,'k-');
    hold on;
    plot(Rv(ind)*rhov,wv(:,ind),'--','color',col);

    xlabel('$r$ [m]');
    ylabel('$w$ [m]');

    figure;
    plot(R*rho,p,'k-');
    hold on;
    plot(Rv(ind)*rhov,pv(:,ind),'--','color',col);
    xlabel('$r$ [m]');
    ylabel('$p$ [Pa]');

end

