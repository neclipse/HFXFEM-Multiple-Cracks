function [ loadrth ] = loadr2r( loadxy, X, Y )
%SIGMAR2X Transform pressure loads form cartesian coordinate system to cylindrical coordinate system,
% i.e., from cartesian traction, Tx and Ty...radial load
[th,~]=cart2pol(X,Y);               % Transform Cartesian coordinates to polar or cylindrical
loadrth=zeros(length(X),2);
for i=1:size(X,1)
dircos=[cos(th(i)),sin(th(i)); -sin(th(i)),cos(th(i))];
T=dircos*loadxy';
loadrth(i,:)=T;
end
end

