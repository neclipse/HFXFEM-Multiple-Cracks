classdef PropagationCheck < handle
properties
    Type            % The type of propagation check
    PVariable       % The principal variable to calculate, either a critical stress or a SIF, etc
    AVariable       % If needed, this is the slot to save the additional value for the check
    Threshold       % Threshold value for the grow, either a critical stress or a SIF, etc
    Tol=0.05;       % Tolerance for the growflag 1<=f<=1+Tol;             
    Growflag=false; % A flag to indicate if the check is passed, false by default
    Mode            % The mode for crack propagation check
end

methods(Abstract)
    cal_pvariable(obj);      % To calculate the principal variable for the crack growth
    growcheck(obj);
end
    
end