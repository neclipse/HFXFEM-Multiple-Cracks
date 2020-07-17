classdef EnrichFun < handle
   properties
       
   end
   
   methods(Abstract)
       enf=calculate(obj);
       enrichelem(obj,elem,id);
       addnodedofs(~,node,id,upindicator)
   end
   
end