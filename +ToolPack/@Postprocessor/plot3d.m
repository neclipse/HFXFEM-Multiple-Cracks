function  h=plot3d( obj, variable, varargin)
%PLOT3D 3D line Plot method of Postprocessor
inc=num2str(obj.Inc);
X=obj.XN;
Y=obj.YN;
var=get(obj,variable);
if ~isempty(varargin)
    var=var(:,varargin{1});                 % To be compatible with EPBar
    component=num2str(varargin{1});
    titlename=['3D plot of ', variable,component,' at nodes after ',inc,' seconds' ];
else
    titlename=['3D plot of ', variable,' at nodes after ',inc,' seconds'];
end
% figure('Name',titlename,'NumberTitle','off');
h=plot3(X,Y,var,'k');
h.LineStyle='none';
h.Marker='.';
title(titlename)
xlabel('x')
ylabel('y')
end





