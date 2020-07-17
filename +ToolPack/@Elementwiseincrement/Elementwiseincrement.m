classdef Elementwiseincrement< ToolPack.CrackGrowIncrementLaw
   properties
       Initial
       Omega   % Proposed new omega for the new crack segment (omega+theta)
       End
   end

   methods
       function obj=Elementwiseincrement(type)
           if nargin>0
               obj.Type=type;
               % single or considering other info to determine if multiple 
               % elements can be inserted at once
           end
       end
       function calincrement(obj)
           % Basically, this is a method of opengeo customized here for
           % flexibility
           % Calculate the end and the increment
           if (strcmp(obj.Type,'single'))
               
               
               
           end
       end
       % can be based on the cohesive zone size and element size
       % should also consider the curvature of the crack
   end
end