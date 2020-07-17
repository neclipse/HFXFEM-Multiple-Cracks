function [psmax,gtheta] = psu( obj )
% 05162019
%PSU principal effective stress update based on obj.Stressp
s=obj.Stressp;     % effective stresses in vector form :[sgmx,sgmy,tauxy,sgmz]
% S=[s(1) s(3) 0;     % effective stress in tensor form
%    s(3) s(2) 0;
%    0   0   s(4)];
%                                                                                   
% in 2D ps1 would be the lastest tensile principal stress
psmax=(s(1)+s(2))/2+sqrt(((s(1)-s(2))/2)^2+s(3)^2);
% ps3=(s(1)+s(2))/2-sqrt(((s(1)-s(2))/2)^2+s(3)^2);
% psmax=max(ps1,ps3);
% the angle between the current global axes to the principal global axes
% if positive then counterclockwisely rotated 
% if negative then clockwisely rotated
% This angle also indicate the rotation angle of crack increment from current crack
gtheta=0.5*atan(2*s(3)/(s(1)-s(2)));    % in radians
% % NO need to do 3D calculation for this 2D plane strain problem
% ps=eig(S);  % the principal stresses are the eigen values of the tensor S;
end

