classdef TractionLaw_bilinearsoftening < EnrichPack.TractionLaw
    % The compression and interpenetration behavior of the cracks are not considered 03102019 
properties
    Stage          % stage of separation of the cohesive crack after solving the equation system
    Stagetrial     % trial stage during crtLHS
    Lambdacr       % Dimensionless parameter to determine if the linear softening begins at the current opening
    AfterPeak=false;      % givepeaktangent is triggered.
end

methods
    function obj = TractionLaw_bilinearsoftening(lambdacr,dcr,tmax,lini,ua,varargin)
        if nargin==4
            super_args={};
        elseif nargin==5
            super_args=ua;
        else
            error('Wrong number of inputs')
        end
        obj=obj@EnrichPack.TractionLaw(super_args);
        obj.Type='bilinear';
        obj.Lambdacr=lambdacr;
        obj.CriticalDisp=dcr;
        obj.PeakTraction=tmax;
%         obj.IniTraction=tini; % not used after 08252019 as initial cohesion is determind from lambdaini
        obj.Lambdaini=lini;
        %% Added on 08152019 to replace the Lambdacr
%         global E;
        agl=ToolPack.GlobalInput();
        E=agl.E;
        if ~isempty(varargin)
            obj.Tnormalc=varargin{1};
            obj.Tshearc=varargin{2};
        else
            obj.Tnormalc=E;
            obj.Tshearc=obj.Mufric*E;
        end
    end
    function stage=checkseparation (obj,mode)
        %Non-dimension shear and normal components
        dshear=obj.Disp(1);
        dnormal=obj.Disp(2);
        dcr=obj.CriticalDisp;
        nds=dshear/dcr;
        ndn=dnormal/dcr;
        if dnormal>=0
            lambdae=sqrt(nds^2+ndn^2);      % nondimensional effective displacement
            % Check the stage of cohesive crack
            if lambdae<=obj.Lambdacr % adhesion stage
                stage=1;
            elseif lambdae<1 % softening stage
                stage=2;
            else % added on 03172019 to ensure the Traction vanish after complete separation.
                % May be modified further to enable healing of kerogen-rich
                % shale.
                stage=3;
                obj.Separated=1;
            end
        else
            lambdae=abs(nds);
            stage=-1;   % contact mode
        end
        % loading or unloading
        if lambdae>=obj.LambdaeL
            % Loading
            obj.LambdaeL=lambdae;
            obj.Loading=1;
        else
            % unloading
            obj.Loading=-1;
        end
        if mode==1           % called when giving tangent for crtLHS
            obj.Stagetrial=stage;
        elseif mode==2
            obj.Stage=stage; % checked after solving equation system, to check if the separation stage match trial
        end
    end

    %% Function declarations
    Tangent=givetangent(obj);
    peakTangent=givepeaktangent(obj);
end
    
end
