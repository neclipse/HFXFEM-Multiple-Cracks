function [ loadxy ] = loadr2x( loadrth, X, Y )
%SIGMAR2X Transform pressure loads form cylindrical coordinate system to cartesian coordinate system,
% i.e., from radial load... to cartesian traction, Tx and Ty...
[th,~]=cart2pol(X,Y);               % Transform Cartesian coordinates to polar or cylindrical
loadxy=zeros(length(X),2);
for i=1:size(X,1)
dircos=[cos(th(i)),-sin(th(i)); sin(th(i)),cos(th(i))];
T=dircos*loadrth';
loadxy(i,:)=T;
end
end

