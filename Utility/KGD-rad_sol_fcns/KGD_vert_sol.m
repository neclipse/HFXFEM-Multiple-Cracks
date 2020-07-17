function [l,w,p,xi,eta]=KGD_vert_sol(Ep,mup,Kp,Cp,Q0,t,N)

xi0=linspace(0,1,N+1);
xi=(xi0(1:end-1)+xi0(2:end))/2;

l=zeros(1,4);
eta=zeros(1,4);
w=zeros(N,4);
p=zeros(N,4);

%M vertex
l(1)=0.6159*(Q0^3*Ep*t^4/mup)^(1/6);
eta(1)=1;
w(:,1)=1.1265*(mup*Q0^3*t^2/Ep)^(1/6)*(1+xi).^0.588.*(1-xi).^(2/3);
p(:,1)=2.7495*(mup*Ep^2/t)^(1/3)*fcn_pres_cal(2/3,0.588,xi',1);
   
%Mt vertex
l(2)=0.3183*Q0*t^(1/2)/Cp;
eta(2)=0;
w(:,2)=0.8165*(mup*Q0^3*t/Ep/Cp^2)^(1/4)*(1+xi).^0.520.*(1-xi).^(5/8);
p(:,2)=3.6783*(Cp^2*mup*Ep^3/t/Q0)^(1/4)*fcn_pres_cal(5/8,0.520,xi',1);
   
%K vertex
l(3)=0.9324*(Ep*Q0*t/Kp)^(2/3);
eta(3)=1;
w(:,3)=0.6828*(Kp^2*Q0*t/Ep^2)^(1/3)*(1-xi.^2).^(1/2);
p(:,3)=0.1831*(Kp^4/Ep/Q0/t)^(1/3)*ones(N,1);
 
%Kt vertex
l(4)=0.3183*Q0*t^(1/2)/Cp;
eta(4)=0;
w(:,4)=0.3989*(Kp^4*Q0^2*t/Ep^4/Cp^2)^(1/4)*(1-xi.^2).^(1/2);
p(:,4)=0.3183*(Kp^4*Cp^2/Q0^2/t)^(1/4)*ones(N,1);

end

