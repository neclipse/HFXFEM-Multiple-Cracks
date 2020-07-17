function xh = fcn_g_del(Kh,Ch)
%Kh - \hat K
%Ch - \hat C

%Kh<0 or complex
iKh=find((Kh<0)|(abs(imag(Kh))>0));
Kh(iKh)=0;
if isempty(iKh)==0
   disp('Warning: Kh is negative or complex in fcn_g_del');
end

%to fix the M vertex
Kh=Kh+eps;

%no propagation in this case
Kh(Kh>1)=1;

%Ch<0 or complex
iCh0=find((Ch<0)|(abs(imag(Ch))>0));
Ch(iCh0)=0;
if isempty(iCh0)==0
   disp('Warning: Ch is negative or complex in fcn_g_del');
end


betam=2^(1/3)*3^(5/6);
betamt=4/15^(1/4)/(sqrt(2)-1)^(1/4);

b0=3*betamt^4/4/betam^3;%b0=0.9912

%function f, solution to differ. equation
f=@(Kh,Ch,C1) (1-Kh.^3-3/2*Ch.*(1-Kh.^2)+3*Ch.^2.*(1-Kh)-3*Ch.^3*2.*atanh((1-Kh)./(2*Ch+1+Kh)) )./(3*C1);

%k-mt edge expanded solution Ch>>1
fkmt=@(Kh,Ch,C1) (1./(4*Ch).*(1-Kh.^4)-1./(5*Ch.^2).*(1-Kh.^5)+1./(6*Ch.^3).*(1-Kh.^6) )./C1;


%functions C1 and C2
C1=@(del) 4*(1-2*del)./(del.*(1-del)).*tan(pi*del);
C2=@(del) 16*(1-3*del)./(3*del.*(2-3*del)).*tan(3*pi/2*del);

%use k-mt edge solution for large values of Ch
iCh=find(Ch>1e3);

del=betam^3/3*f(Kh,b0*Ch,betam^3/3).*(1+b0*Ch);
del(iCh)=betam^3/3*fkmt(Kh(iCh),b0*Ch(iCh),betam^3/3).*(1+b0*Ch(iCh));

del(del<=0)=1e-6;
del(del>=1/3)=1/3-1e-6;

bh=C2(del)./C1(del);

%delta-correction
xh=f(Kh,Ch.*bh,C1(del));
xh(iCh)=fkmt(Kh(iCh),Ch(iCh).*bh(iCh),C1(del(iCh)));

end

