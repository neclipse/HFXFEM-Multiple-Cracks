function discretize(obj,np)                                 % discretize the smooth curve with the given function handle and limits
ind=1:np;
description=obj.Description;
if obj.Linetype == 2
    % with a function handle
    Funhandle=description{1};
    limits=description{3};
    % cut the line into pieces
    indpoints=linspace(limits(1),limits(2),np);                      % independent variable x or y
    dpoints=Funhandle(indpoints);                                % dependent variable y or x
    if description{2} == 1                                           % The independent variable for limits is x
        obj.Segments=transpose([ind;indpoints;dpoints]);             % Segments [ind,x,y]
    elseif description{2} == 2                                       % The independent variable for limits is y
        obj.Segments=transpose([ind;dpoints;indpoints]);             % Segments [ind,x,y]
    end
elseif obj.Linetype ==1
    % with discrete point list
    [xlist,I]=sort(obj.Segments(:,2));
    obj.Segments=obj.Segments(I,:);
    % As the tip direction, the sign of phi in opengeo.initiate and ntaud
    % of line in opengeo.intersection all assumes that the order of x is
    % ascending, we had better enforce this condition to guarantee the
    % consistency.
    ylist=obj.Segments(:,3);
    if all(xlist-xlist(1)==0)
        limits=[obj.Segments(1,3),obj.Segments(end,3)];               % ylimits of the segments
        indpoints=linspace(limits(1),limits(2),np);                   % independent variable y
        dpoints=interp1(ylist,xlist,indpoints);
        obj.Segments=transpose([ind;dpoints;indpoints]);
    else
        limits=[obj.Segments(1,2),obj.Segments(end,2)];               % xlimits of the segments
        indpoints=linspace(limits(1),limits(2),np);                   % independent variable y
        dpoints=interp1(xlist,ylist,indpoints);
        obj.Segments=transpose([ind;indpoints;dpoints]);
    end
end
end

