function h=plotsurf( obj, variable, meshstructure,deformedflag,varargin )
%PLOTCONTOUR Method of Postprocessor class
inc=num2str(obj.Inc);
var=get(obj,variable);
if ~isempty(varargin)
    var=var(:,varargin{1});                 % To be compatible with EPBar
    component=num2str(varargin{1});
    titlename=['Surface plot of the ', variable,component,' averaged at nodes after ',inc,' seconds'];
else
    titlename=['Surface plot of the ', variable,' averaged at nodes after ',inc,' seconds'];
end
%reshape the vectors to be compatible with the requirement of the built-in 'surf'
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
h=surf(X,Y,var);
% daspect([1000 1000 1])
h.EdgeColor='none';
colormap (jet)
colorbar;
% title(titlename)
% set(get(gca,'title'),'Position',[0.5,5,1]);
xlabel('x')
ylabel('y')
% view(0,90)
end
