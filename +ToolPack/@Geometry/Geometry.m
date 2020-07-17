classdef Geometry < handle
   properties
       Geotype                                                              % 1. (open end) piece-linear curve 2. (closed end) analytical function like a circle, an ellipse,etc
       Mesh
       Nodedict
       Elemdict
   end
   properties (Access = protected)
       
   end
   methods(Abstract)
       initiate(obj)                                                        % to find the inititial interacted elements with the geometry 
       update(obj)                                                          % update the geometry and the interacted elements with the updated geometry
       interactwith(obj)                                                    
   end
   methods
       function obj=Geometry(mesh,nodedict,elemdict)
          if nargin>0
          obj.Mesh=mesh;
          obj.Nodedict=nodedict;
          obj.Elemdict=elemdict;
          end
       end
       elems=findelems(obj,plist,varargin);                        % find the interacting elements with a list of plist(x,y;...) with options 'inside', 'edge', 'in_edge'
   end

end