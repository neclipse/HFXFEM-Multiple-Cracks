classdef TractionLaw_unified < EnrichPack.TractionLaw
    % The compression and interpenetration behavior of the cracks are not considered 03102019 
properties
    Stage          % stage of separation of the cohesive crack after solving the equation system
    Stagetrial     % trial stage during crtLHS
    Lambdacr       % Dimensionless parameter to determine if the linear softening begins at the current opening
    KrgTraction    % The turning value corresponding to Lambdacr
    AfterPeak=false;      % givepeaktangent is triggered.
end

methods
    function obj = TractionLaw_unified(lambdacr,dcr,tkrg,tini,lini,gbinp)
        super_args={};
        obj=obj@EnrichPack.TractionLaw(super_args);
        obj.Type='unified';
        obj.Lambdacr=lambdacr;
        obj.CriticalDisp=dcr;
        obj.KrgTraction=tkrg;
        obj.IniTraction=tini;
        obj.PeakTraction=max(tini,tkrg);
        obj.Lambdaini=lini;
%        obj.Lambdaini=lambdaini;% not used FOR THIS CLASS AS LAMBDAINI
%        SHOULD BE ZERO.
        %% Added on 08152019 to replace the Lambdacr
        obj.Tnormalc=gbinp.E*1e2; % * 1e3 Increased to a very big number to reduce penetration 01/2021.
        obj.Tshearc=obj.Mufric*gbinp.E;
    end
    function stage=checkseparation (obj,mode)
        % This function remain unused for the current scheme 052620
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
