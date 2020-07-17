function Tangent=givetangent(obj)
%Non-dimension shear and normal components
dshear=obj.Disp(1);
% In case of Mode I crack, shear displacement is minmal. It is
% better to drop tiny shear displacement and related shear traction
% which are mainly caused by numerical error. 07072019
if abs(dshear)<1e-12
    dshear=0;
end
dnormal=obj.Disp(2);
lcr=obj.Lambdacr;
dcr=obj.CriticalDisp;
scr=obj.PeakTraction;
nds=dshear/dcr;
ndn=dnormal/dcr;
if dnormal>=0
    % The tensile mode
    lambdae=sqrt(nds^2+ndn^2);      % nondimensional effective displacement
    % Assume zero crack stiffness if crack opening is zero, however, a
    % sudden jump in crack stiffness will be observed when lamdae
    % increases from zero to a tiny value. NEED TO BE DISCUSSED LATER
    % 02/04/2019. Maybe it is the reason why Khoei in 2010 used
    % bilinear cohesive law is to avoid such discontinuity.
    if lambdae<=obj.Lambdacr % adhesion stage
        if lambdae>=obj.LambdaeL
            obj.LambdaeL=lambdae;
            obj.Loading=1;
        else
            obj.Loading=-1;
        end
        % loading and unloading is the same
        Tss=scr/lcr/dcr;
        Tsn=0;
        Tnn=Tss;
        Tangent=[Tss,Tsn;Tsn,Tnn];
        obj.Tangent=Tangent;
    elseif lambdae<1 % softening stage
        if lambdae>=obj.LambdaeL 
            % Loading or just afterpeak
            obj.LambdaeL=lambdae;
            obj.Loading=1;
            temp1=dcr*scr/(1-lcr);
            Tss=-temp1*(dshear/lambdae/dcr^2)^2+...
                (1-lambdae)*temp1*(1/lambdae/dcr^2-dshear^2/lambdae^3/dcr^4);
            Tsn=-temp1/lambdae^3*dshear*dnormal/dcr^4;
            Tnn=-temp1*(dnormal/lambdae/dcr^2)^2+...
                (1-lambdae)*temp1*(1/lambdae/dcr^2-dnormal^2/lambdae^3/dcr^4);
        else
            % unloading
            obj.Loading=-1;
            Tss=scr/dcr*(1-obj.LambdaeL)/(1-lcr)/obj.LambdaeL;
            Tsn=0;
            Tnn=Tss;
        end
    else % added on 03172019 to ensure the Traction vanish after complete separation.
        % May be modified further to enable healing of kerogen-rich
        % shale.
        obj.Separated=1;
        Tss=0; Tsn=0; Tnn=0;
    end
else
    %0825 Add the compressive mode, (temporary relationship)
    lambdae=abs(nds);      % nondimensional effective displacement has only shear component
    if lambdae>=obj.LambdaeL
        % Loading
        obj.LambdaeL=lambdae;
        obj.Loading=1;
    else
        % unloading
        obj.Loading=-1;
    end
    Tss=obj.Tshearc; Tsn=0; Tnn=obj.Tnormalc; 
end
Tangent=[Tss,Tsn;Tsn,Tnn];
obj.Tangent=Tangent;
end

