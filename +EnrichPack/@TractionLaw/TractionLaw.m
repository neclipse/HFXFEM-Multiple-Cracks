classdef TractionLaw < handle
properties
    Type
    Disp         % Relative displacement or displacement discontinuity in local orthogonal coordinate system
    Loading=1;   % Flag to indicate loading(1) or unloading(-1)
    Tangent      % Tanget matrix in terms of local orthogonal coordinate system
    LambdaeL=0;    % The maximum effetive separation in history
    Separated    % A flag to know if the crack surfaces have been completely separated.
    IniTraction=0;  % Initial traction at the initial state (may not be zero opening)
    CriticalDisp % Critical displacement where complete separation is in place (m) 
    PeakTraction % Critical displacement where Tensile(Positive) traction reaches its peak (Mpa) 
    Tnormalc        % normal modulus for compressive contact mode, approximately is Young's modulus although it is not strain based law
    Tshearc         % shear modulus for compressive contact mode, approximately is internal friction* Young's modulus before breaking Mohr-Coulomb law.
    Lambdaini=0.05; % 0.05
    Mufric=0.3;          % internal friction coefficient in Mohr-Coulomb law.
    Slipflag=false;  % if |tsc|+Mufric*|tnc| is not met
end
methods
    function obj=TractionLaw(varargin)
        if nargin>0

        end
    end
end
methods(Abstract)
    givetangent(obj)
end
    
end
