function Deltap = fcn_Delta_p(Kh,Ch,p)
%Kh - \hat K
%Ch - \hat C
%p - parameter from 0 to 1

%Kh<0 or complex
iKh=find((Kh<0)|(abs(imag(Kh))>0));
Kh(iKh)=0;
if isempty(iKh)==0
   disp('Warning: Kh is negative or complex in fcn_Delta_p');
end

%to fix the M vertex
Kh=Kh+eps;

%no propagation in this case
Kh(Kh>1)=1;

%Ch<0 or complex
iCh0=find((Ch<0)|(abs(imag(Ch))>0));
Ch(iCh0)=0;
if isempty(iCh0)==0
   disp('Warning: Ch is negative or complex in fcn_Delta_p');
end


betam=2^(1/3)*3^(5/6);
betamt=4/15^(1/4)/(sqrt(2)-1)^(1/4);

b0=3*betamt^4/4/betam^3;%b0=0.9912

%function f, solution to diffier. equation
f=@(Kh,Ch,C1) (1-Kh.^3-3/2*Ch.*(1-Kh.^2)+3*Ch.^2.*(1-Kh)-3*Ch.^3*2.*atanh((1-Kh)./(2*Ch+1+Kh)) )./(3*C1);

%k-mt edge expanded solution
fkmt=@(Kh,Ch,C1) (1./(4*Ch).*(1-Kh.^4)-1./(5*Ch.^2).*(1-Kh.^5)+1./(6*Ch.^3).*(1-Kh.^6) )./C1;

%use k-mt edge solution for large values Ch
iCh=find(Ch>1e3);

Delta=betam^3/3*f(Kh,Ch*b0,betam^3/3).*(1+b0*Ch);
Delta(iCh)=betam^3/3*fkmt(Kh(iCh),Ch(iCh)*b0,betam^3/3).*(1+b0*Ch(iCh));

Deltap=(1-p+p*f(Kh,Ch*b0,betam^3/3).*(betam^3+betamt^4*Ch)).*Delta;
Deltap(iCh)=(1-p+p*fkmt(Kh(iCh),Ch(iCh)*b0,betam^3/3).*(betam^3+betamt^4*Ch(iCh))).*Delta(iCh);

end

