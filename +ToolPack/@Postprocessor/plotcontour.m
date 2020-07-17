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
var=transpose(reshape(var,meshstructure));
if deformedflag     % surfplot on the deformed mesh
    dispmin=min(min([obj.UX,obj.UY]));
    dispmax=max(max([obj.UX,obj.UY]));
    maxabsdisp=max(abs(dispmin),abs(dispmax));
    scale=10^(abs(fix(log10(maxabsdisp)))+0.5);      % scaling factor
    Xd=obj.XN+obj.UX*scale;
    Yd=obj.YN+obj.UY*scale;
    X=transpose(reshape(Xd,meshstructure));
    Y=transpose(reshape(Yd,meshstructure));
else                % surfplot on the undeformed mesh
    X=transpose(reshape(obj.XN,meshstructure));
    Y=transpose(reshape(obj.YN,meshstructure));
end
figure('Name',titlename,'NumberTitle','off');
h=contourf(X,Y,var,ncontour);
daspect([1 1 1])
colormap (jet)
colorbar
% title(titlename);
xlabel('x');
ylabel('y');
zlabel(variable)
end

