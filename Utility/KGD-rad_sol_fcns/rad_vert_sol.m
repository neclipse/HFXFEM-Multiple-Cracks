function [R,w,p,rho,eta]=rad_vert_sol(Ep,mup,Kp,Cp,Q0,t,N)

rho0=linspace(0,1,N+1);
rho=(rho0(1:end-1)+rho0(2:end))/2;
    
R=zeros(1,4);
eta=zeros(1,4);
w=zeros(N,4);
p=zeros(N,4);

%M vertex
R(1)=0.6944*(Q0^3*Ep*t^4/mup)^(1/9);
eta(1)=1;
w(:,1)=1.1901*(mup^2*Q0^3*t/Ep^2)^(1/9)*(1+rho).^0.487.*(1-rho).^(2/3);
p(:,1)=2.4019*(mup*Ep^2/t)^(1/3)*fcn_pres_cal(2/3,0.487,rho',2);
   
%Mt vertex
R(2)=0.4502*(Q0^2*t/Cp^2)^(1/4);
eta(2)=0;
w(:,2)=1.0574*(mup^4*Q0^6*t/Ep^4/Cp^2)^(1/16)*(1+rho).^0.397.*(1-rho).^(5/8);
p(:,2)=3.0931*(Cp^6*mup^4*Ep^(12)/t^3/Q0^2)^(1/16)*fcn_pres_cal(5/8,0.397,rho',2);
   
%K vertex
R(3)=0.8546*(Ep*Q0*t/Kp)^(2/5);
eta(3)=1;
w(:,3)=0.6537*(Kp^4*Q0*t/Ep^4)^(1/5)*(1-rho.^2).^(1/2);
p(:,3)=0.3004*(Kp^6/Ep/Q0/t)^(1/5)*ones(N,1);
 
%Kt vertex
R(4)=0.4502*(Q0^2*t/Cp^2)^(1/4);
eta(4)=0;
w(:,4)=0.4744*(Kp^8*Q0^2*t/Ep^8/Cp^2)^(1/8)*(1-rho.^2).^(1/2);
p(:,4)=0.4139*(Kp^8*Cp^2/Q0^2/t)^(1/8)*ones(N,1);

end

