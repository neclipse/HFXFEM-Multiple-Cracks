function [xyl,r]=callocal(obj,x,y)                                   % calculate the local coordinates (xl,yl) of a point (x,y), and the distance from the tip
% make sure x and y are row vectors.
if size(x,1)>1
    x=x';
end
if size(y,1)>1
    y=y';
end
o=obj.Omega;
T=[cos(o) sin(o); -sin(o) cos(o)];                               % transformation matrix
xc=x-obj.Xct;
yc=y-obj.Yct;
xyl=T*[xc;yc];
r=sqrt(xyl(1,:).^2+xyl(2,:).^2);
end