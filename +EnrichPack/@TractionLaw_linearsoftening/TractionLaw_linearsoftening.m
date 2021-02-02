classdef TractionLaw_linearsoftening < EnrichPack.TractionLaw
    properties
        AfterPeak=false;
    end
    methods
        function obj = TractionLaw_linearsoftening(dcr,tmax,lambdaini,gbinp)
            super_args={};
            obj=obj@EnrichPack.TractionLaw(super_args);
            obj.Type='linear';
            obj.CriticalDisp=dcr;
            obj.PeakTraction=tmax;
%             obj.IniTraction=tini; % not assigned after 08252019 as initial cohesion is determind from lambdaini
            obj.Lambdaini=lambdaini;
            % Added on 08152019 to replace the Lambdacr
            obj.Tnormalc=gbinp.E*1e3;
            obj.Tshearc=obj.Mufric*gbinp.E;
        end
        peaktangent=givepeaktangent(obj);
        Tangent=givetangent(obj);
          
    end
    
end
