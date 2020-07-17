function [h,S,var]= plotxy( obj,variable,indices,type,r0,varargin)
%PLOTXY X-Y Plot method of Postprocessor
%   To plot variable versus distant at nodes/element(one point per element) defined by the indices 
inc=num2str(obj.Inc);
X=obj.XN(indices);
Y=obj.YN(indices);
S=sqrt((X-X(1)).^2+(Y-Y(1)).^2);                % Distance along the path measured from the origin
var=get(obj,variable);
if nargin > 5
    var=var(:,varargin{1});                 % To be compatible with EPBar
    component=num2str(varargin{1});
    titlename=['X-Y plot of ', variable,component,' averaged at nodes after ',inc,' seconds'];
else
    titlename=['X-Y plot of ', variable,' averaged at nodes after ',inc,' seconds'];
end
S=S/r0;                                     % dimensionless distance
var=var(indices);
figure('Name',titlename,'NumberTitle','off');
switch type
    case 'linear'
        h=plot(S,var,'k');
    case 'xlog'
        h=semilogx(S,var,'k');
    case 'ylog'
        h=semilogy(S,var,'k');
end
% h.Marker='+';
% title(titlename)
% xlabel('Dimensionless distance along the path')
% ylabel(variable)
end


