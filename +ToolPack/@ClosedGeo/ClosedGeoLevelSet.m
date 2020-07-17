classdef ClosedGeoLevelSet < ToolPack.Geometry
   properties (Access = public)
       Nodedict
       Funhandle
   end
   properties (Access = protected)
       
   end
   methods
       initiate(obj)                                                        % to find the inititial interacted elements with the geometry 
       update(obj)                                                          % update the geometry and the interacted elements with the updated geometry
       interactwith(obj,elem)
   end
   methods
       
   end
end