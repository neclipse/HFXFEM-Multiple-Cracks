curvexy = [0 0;1 0;2 1;3 2];
mapxy = [3 4;.5 .5;3 -1;4 1];
[xy,distance,t,minindices] = distance2curve(curvexy,mapxy,'linear');
% xy =
%                          2                         1
%          0.470588235294118         0.617647058823529
%                        1.5                       0.5
% distance =
%           3.16227766016838
%          0.121267812518166
%           2.12132034355964
% t =
%          0.485194315877587
%          0.802026225550702
%           0.34308419095021


plot(curvexy(:,1),curvexy(:,2),'k-o',mapxy(:,1),mapxy(:,2),'r*')
hold on
plot(xy(:,1),xy(:,2),'g*')
line([mapxy(:,1),xy(:,1)]',[mapxy(:,2),xy(:,2)]','color',[0 0 1])
axis equal