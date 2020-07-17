function painter( obj,varargin )
%PAINTER Method of class 'postprocessor' to plot the results
%   This method can provide various visualizationo of the results calculted
%   by the NewtonRaphson method
%% --   Syntax:
%   painter(obj)                        Default painter: plot 3d figures of all common variables at once
%   painter(obj,meshstructure)          Advanced painter: plot 2d filled contour figures for all common variables when mesh is structured
%   painter(obj,mode,indices)           User defined X-Y plot: plot X-Y figures of all variables along specified node path or element path according to the mode
%% Painting according to user's request
if nargin==1                            % Default paint: plot all suface figures at once
    obj.plot3d('UX');
    obj.plot3d('UY');
    obj.plot3d('EPBar');
    obj.plot3d('SG',1);
    obj.plot3d('SG',2);
elseif nargin==3
    % Mesh structure: idim*jdim which is exclusively for structured mesh, the nodes/elements are continuously ordered along a j-direction
    % if nodes/elements are continuously ordered along i-direction, then meshstructure is jdim*idim
    nodestructure=varargin{1}; 
    pointstructure=varargin{2};
    obj.plotcontour('UX',nodestructure,0);
    obj.plotcontour('UY',nodestructure,0);
    obj.plotcontour('EPBar',pointstructure,0);
    obj.plotcontour('SG',pointstructure,0,1);
    obj.plotcontour('SG',pointstructure,0,2);
elseif nargin==4
    nodestructure=varargin{1}; 
    pointstructure=varargin{2};
    nocontour=varargin{3};
    obj.plotcontour('UX',nodestructure,nocontour);
    obj.plotcontour('UY',nodestructure,nocontour);
    obj.plotcontour('EPBar',pointstructure,1,nocontour);
    obj.plotcontour('SG',pointstructure,0,1,nocontour);
    obj.plotcontour('SG',pointstructure,0,2,nocontour);
    
end
end

