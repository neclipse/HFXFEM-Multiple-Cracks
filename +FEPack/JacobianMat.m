classdef JacobianMat < handle
    properties
        Id
        Locarray
        LocarrayEnr
        LocarrayU
        LocarrayUEnr
        LocarrayP
        LocarrayPEnr
        LocarrayAll
        IntLoadVecAll             % Internal load vector corresponding to the standard and the current enrichment item
        DUn2Enr                    % incremental change of Uenr, dUn2i+dUn2(i-1)+...
        Un2iEnr                    % Un1+dUn2i
        Pn2iEnr                    % Pn1+dPu2i
        Un2t1Enr
        Un2t2Enr
        Pn2t1Enr
		F_coh_old                 % Converged Cohesion force vector averaged to the nodes at last increment==zeros(8,1)
		F_coh_i					   % Trial Cohesion force vector averaged to the nodes at the current iteration i
        Musus
        Musue
        Mueus
        Mueue
        Kusus
        Kusue
        Kueus
        Kueue
        Qusps
        Quspe
        Queps
        Quepe
        Hpsps
        Hpspe
        Hpeps
        Hpepe
        Spsps
        Spspe
        Speps
        Spepe
        WJ                          % The whole Jacobian Matrix for the current element;
        WJnew                       % The whole Jacobian Matrix only used when it moves to the next increment.
		Kc 							% Cohesion Tangent matrix in global coordinate system
		Qintueps
		Qintuepe
		Hintpsps
		Hintpspe
		Hintpeps
		Hintpepe
		Sintpsps
		Sintpspe
		Sintpeps
		Sintpepe
    end
    
    methods
        function obj=JacobianMat(Id)
            if nargin>0
                obj.Id=Id;
            end
            
        end
    end
end