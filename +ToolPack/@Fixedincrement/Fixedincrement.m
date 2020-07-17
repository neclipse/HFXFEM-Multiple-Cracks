classdef Fixedincrement< ToolPack.CrackGrowIncrementLaw
   properties
       
   end

   methods
       function obj=Fixedincrement(type)
           if nargin>0
               obj.Type=type;
           end
       end
       function calincrement(obj,x,varargin)
           obj.Increment=x;
       end
       % can be based on the cohesive zone size and element size
       % should also consider the curvature of the crack
   end
end