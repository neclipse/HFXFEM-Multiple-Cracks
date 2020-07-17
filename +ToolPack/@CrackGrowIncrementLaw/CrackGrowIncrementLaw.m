classdef CrackGrowIncrementLaw < handle
   properties
       Increment        % The crack grow increment length
       Type             % Paris, User, Elementwise, etc
       Elemahead
   end
   
   
   methods(Abstract)
       calincrement(obj);
       % can be based on the cohesive zone size and element size
       % should also consider the curvature of the crack
   end
%    methods
% 
%    end
    
    
end