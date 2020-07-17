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
scr=obj.KrgTraction;
sini=obj.IniTraction;
nds=dshear/dcr;
ndn=dnormal/dcr;
kfirst=(scr-sini)/lcr;
if dnormal>=0
    % The tensile mode
    lambdae=sqrt(nds^2+ndn^2);      % nondimensional effective displacement
    if obj.Lambdacr>0 && lambdae<=obj.Lambdacr % first stage
        obj.Stage=1;
%         if lambdae>=obj.LambdaeL*0.98
            obj.LambdaeL=lambdae;
            obj.Loading=1;
            % Trick to avoid infinity for the initial zero separation case
            if dnormal==0 && dshear==0
                Tss=kfirst/dcr;
                Tsn=0;
                Tnn=kfirst/dcr;
            else
                Tss=kfirst/dcr+dnormal^2*sini/dcr^3/lambdae^3;
                Tsn=-dnormal*dshear*sini/(dcr*lambdae)^3;
                Tnn=kfirst/dcr+dshear^2*sini/dcr^3/lambdae^3;
            end
%         else
%             obj.Loading=-1;
% 			Tss=(khard+sini/obj.LambdaeL)/dcr;
% 			Tsn=0;
% 			Tnn=Tss;
%         end
    elseif lambdae<1 % second stage
        obj.Stage=2;
%         if lambdae>=obj.LambdaeL*0.98
            % Loading or just afterpeak
            obj.LambdaeL=lambdae;
            obj.Loading=1;
            temp1=dcr*scr/(1-lcr);
            % Trick to avoid infinity for the initial zero separation case
            if dnormal==0 && dshear==0
                Tss=-scr/dcr;
                Tsn=0;
                Tnn=-scr/dcr;
            else
                Tss=-temp1*(dshear/lambdae/dcr^2)^2+...
                    (1-lambdae)*temp1*(1/lambdae/dcr^2-dshear^2/lambdae^3/dcr^4);
                Tsn=-temp1/lambdae^3*dshear*dnormal/dcr^4;
                Tnn=-temp1*(dnormal/lambdae/dcr^2)^2+...
                    (1-lambdae)*temp1*(1/lambdae/dcr^2-dnormal^2/lambdae^3/dcr^4);
            end
%         else
%             % unloading
%             obj.Loading=-1;
%             Tss=scr/dcr*(1-obj.LambdaeL)/(1-lcr)/obj.LambdaeL;
%             Tsn=0;
%             Tnn=Tss;
%         end
    else % added on 03172019 to ensure the Traction vanish after complete separation.
        % May be modified further to enable healing of kerogen-rich
        % shale.
        obj.Stage=3;
        obj.Separated=1;
        Tss=0; Tsn=0; Tnn=0;
    end
else
    %0825 Add the compressive mode, (temporary relationship)
    lambdae=abs(nds);      % nondimensional effective displacement has only shear component
    obj.Stage=-1;
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

