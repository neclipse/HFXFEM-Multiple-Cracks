function Tangent=givetangent(obj)
%Non-dimension shear and normal components
dshear=obj.Disp(1);
dnormal=obj.Disp(2);
dcr=obj.CriticalDisp;
scr=obj.PeakTraction;
nds=dshear/dcr;
ndn=dnormal/dcr;
% Assume zero crack stiffness if crack opening is zero, however, a
% sudden jump in crack stiffness will be observed when lamdae
% increases from zero to a tiny value. NEED TO BE DISCUSSED LATER
% 02/04/2019. Maybe it is the reason why Khoei in 2010 used
% bilinear cohesive law is to avoid such discontinuity.
if dnormal>=0
    % The tensile mode
    lambdae=sqrt(nds^2+ndn^2);      % nondimensional effective displacement
    if lambdae<0.05
        % initially compressive
        obj.LambdaeL=lambdae;
        obj.Loading=1;
%         Tss=-scr/dcr;        % an arbitray large positive number, changed to zero
%         Tsn=0;
%         Tnn=scr/dcr/obj.Lambdacr;   % initialize for newly generated cohesive segments
        % initially tensile
        Tss=-scr/dcr;
        Tsn=0;
        Tnn=-scr/dcr;
    elseif lambdae<1
%         if lambdae>=obj.LambdaeL*0.99 % to avoid numerical unstability induced unloading.
            % Loading
            obj.LambdaeL=lambdae;
            obj.Loading=1;
            % Trick to avoid infinity for the initial zero separation case
            if dnormal==0 && dshear==0
                % the limit value when dnormal=0
                Tss=-scr/dcr;
                Tsn=0;
                % the limit value when dshear=0
                Tnn=-scr/dcr;
            else
                Tss=-dshear^2*scr/dcr^3/lambdae^2+...
                    (1-lambdae)*dnormal^2*scr/dcr^3/lambdae^3;
                Tsn=-scr*dshear*dnormal/dcr^3/lambdae^3;
                Tnn=-dnormal^2*scr/dcr^3/lambdae^2+...
                    (1-lambdae)*dshear^2*scr/dcr^3/lambdae^3;
            end
%         else 
%             unloading
%             obj.Loading=-1;
%             Tss=scr/dcr*(1-obj.LambdaeL)/obj.LambdaeL;
%             Tsn=0;
%             Tnn=Tss;
%         end
    else
        % added on 07042019 to ensure the Traction vanish after complete separation.
        % May be modified further to enable healing of kerogen-rich
        % shale.
        obj.Separated=1;
        Tss=0; Tsn=0; Tnn=0;
    end
else
    %0728 Add the compressive mode, (temporary relationship)
    lambdae=abs(nds);      % nondimensional effective displacement has only shear component
%     if lambdae>=obj.LambdaeL
%         % Loading
%         obj.LambdaeL=lambdae;
%         obj.Loading=1;
%     else
%         % unloading
%         obj.Loading=-1;
%     end
%     Tss=obj.Tshearc; Tsn=0; Tnn=obj.Tnormalc;
    %% The tractionlaw in compression contact mode based on separation is
    % deemed not realistic and abandoned on 08/15/2019
    if lambdae<1
%         if lambdae>=obj.LambdaeL
            % Loading
            obj.LambdaeL=lambdae;
%             obj.Loading=1;
            Tss=-scr/dcr;
            Tsn=0;
            Tnn=obj.Tnormalc; % to prevent interpenetration
%         else
%             % unloading
%             obj.Loading=-1;
%             Tss=scr/dcr*(1-obj.LambdaeL)/obj.LambdaeL;
%             Tsn=0;
%             Tnn=obj.Tnormalc;
%         end
    else
        % added on 07042019 to ensure the Traction vanish after complete separation.
        % May be modified further to enable healing of kerogen-rich
        % shale.
        obj.Separated=1;
        Tss=obj.Tshearc; Tsn=0; Tnn=obj.Tnormalc;
    end
end
% For near K mode, only focus on mode 1.
Tangent=[Tss,Tsn; Tsn,Tnn];
obj.Tangent=Tangent;
end




