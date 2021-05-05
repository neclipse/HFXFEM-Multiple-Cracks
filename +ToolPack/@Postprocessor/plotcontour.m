function h=plotcontour( obj, variable, meshstructure,deformedflag,varargin )
%PLOTCONTOUR Method of Postprocessor class
ncontour=10;                                  % default number of contour lines
if nargin==5
    component=varargin{1};
elseif nargin==6
    component=varargin{1};
    ncontour=varargin{2};                    % specify the number of contour lines
end
inc=num2str(obj.Inc);
var=get(obj,variable);
if nargin>4
    var=var(:,component);                 % To be compatible with EPBar
    component=num2str(component);
    titlename=['Countour plot of the ', variable,component,' averaged at nodes after ',inc,' seconds'];
else
    titlename=['Countour plot of the ', variable,' averaged at nodes after ',inc,' seconds'];
end
if deformedflag     % surfplot on the deformed mesh
    dispmin=min(min([obj.UX,obj.UY]));
    dispmax=max(max([obj.UX,obj.UY]));
    maxabsdisp=max(abs(dispmin),abs(dispmax));
    scale=10^(abs(fix(log10(maxabsdisp)))+0.5);      % scaling factor
    Xd=obj.XN+obj.UX*scale;
    Yd=obj.YN+obj.UY*scale;
else                % surfplot on the undeformed mesh
    Xd=obj.XN;
    Yd=obj.YN;
end
figure('Name',titlename,'NumberTitle','off');
if length(meshstructure)==2
    var=transpose(reshape(var,meshstructure));
    X=transpose(reshape(Xd,meshstructure));
    Y=transpose(reshape(Yd,meshstructure));
    h=contourf(X,Y,var,ncontour);
else 
    tri = delaunay(Xd,Yd);
    h = tricontour(tri, Xd, Yd, var*1000, ncontour);
end
daspect([1 1 1])
colormap (jet)
colorbar
% title(titlename);
xlabel('x');
ylabel('y');
zl=strcat(variable," (MPa)");
zlabel(zl)
end

