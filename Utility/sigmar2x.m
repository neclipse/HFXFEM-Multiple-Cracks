function [ sigmax ] = sigmar2x( sigmar, X, Y )
%SIGMAR2X Transform stress form cylindrical coordinate system to cartesian coordinate system,
% i.e., from radial and tangent stress... to sigmax and sigmay...
%   For the plane strain plastic problem, the stress array has four
%   components.
[th,~]=cart2pol(X,Y);                                   % Transform Cartesian coordinates to polar or cylindrical
sigmarth=sigmar(1:3);
sigmax=zeros(length(X),4);
sigmax(:,4)=sigmar(4);
for i=1:size(X,1)
T=[cos(th(i))^2, sin(th(i))^2,-sin(2*th(i));
    sin(th(i))^2,cos(th(i))^2,sin(2*th(i));
    0.5*sin(2*th(i)),-0.5*sin(2*th(i)),cos(2*th(i))];
sigma=T*sigmarth';
sigmax(i,1:3)=sigma';
end
end

