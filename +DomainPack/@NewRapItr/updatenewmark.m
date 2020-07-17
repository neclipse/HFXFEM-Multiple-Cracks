function updatenewmark(obj,varargin)
if ~isempty(varargin)
    if nargin==2
        mode=1;% GN11
        theta=varargin{1};
        
    elseif nargin==4
        mode=2;% GN22 for u and GN 11 for p
        gamma=varargin{1};
        beta=varargin{2};
        theta=varargin{3};
    else
        error('wrong input numbers for the newmark scheme');
    end
else
    % The selection of Newmark parameters requires attention (03262019)
    mode=2;
    gamma=0.5;   % 1/2-> second order accurate.
    beta=0.25;   % 1/4-> constant acceleration; 1/6_> linear acceleration
    theta=0.75;  % was 0.5 and the oscillation in pressure was too strong. changed to 1 on 03262019
end
if obj.IInc>1
    dt=(obj.Timeinc(obj.IInc)-obj.Timeinc(obj.IInc-1));
else
    dt=obj.Timeinc(1);
end
obj.Dt=dt;
switch mode
    case 1 % GN11
        a1p=1/(theta*dt);
        a3p=1/theta-1;
        a0=0; a2=0; a4=0;
        a1=a1p; a3=a3p;a5=0;
    case 2   % GN22 and GN11
        a0=1/(beta*dt^2);
        a1=gamma/(beta*dt);
        a2=1/(beta*dt);
        a3=gamma/beta-1;
        a4=1/(2*beta)-1;
        a5=dt*(gamma/(2*beta)-1);
        a1p=1/(theta*dt);
        a3p=1/theta-1;
end
obj.Newmark=struct('a0',a0,'a1',a1,'a2',a2,'a3',a3,'a4',a4,'a5',a5,'a1p',a1p,'a3p',a3p,'dt',dt);
end
