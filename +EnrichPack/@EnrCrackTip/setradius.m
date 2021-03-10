function setradius(obj,varargin)
if isempty(varargin)
    % By default, the interaction radius is three times of the
    % typical element size in the mesh around the crack tip.
    ratio=3;
else
    ratio=varargin{1};
end
baselength=max(obj.INTELEM.Length);
length= ratio*baselength;
obj.Radius=length;
end

