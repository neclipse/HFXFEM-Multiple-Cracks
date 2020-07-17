classdef CrackGrowDirectionLaw < handle
   properties
       Theta=0;        % counterclockwise deflection angle from the current crack tip for the next step
       Omega        % current tip angle from x-axis in global coordinate system
       PVariable    % can be effective stress vector or stress intensity factor
   end
   
   
   methods(Abstract)
       theta = growdirection(obj);
   end
    
    
end