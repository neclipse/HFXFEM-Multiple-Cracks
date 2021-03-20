function [psmax,gtheta] = psu_exact(  obj, s,  tip, omega )
% 05162019
%PSU principal effective stress update based on obj.Stressp
% effective stresses in vector form :[sgmx;sgmy;tauxy;sgmz]
% one column stands for one entry      
% S=[s(1) s(3) 0;     % effective stress in tensor form
%    s(3) s(2) 0;
%    0   0   s(4)];
%                                                                                   
% % in 2D ps1 would be the lastest tensile principal stress
% % psmax=(s(1,:)+s(2,:))/2+sqrt(((s(1,:)-s(2,:))/2).^2+s(3,:).^2);
% % ps3=(s(1)+s(2))/2-sqrt(((s(1)-s(2))/2)^2+s(3)^2);
% % psmax=max(ps1,ps3);
% % the angle between the current global axes to the principal global axes
% % if positive then counterclockwisely rotated 
% % if negative then clockwisely rotated
% % This angle also indicate the rotation angle of crack increment from current crack
% % gtheta=0.5*atan(2*s(3,:)./(s(1,:)-s(2,:)));    % in radians
% % % NO need to do 3D calculation for this 2D plane strain problem
% % ps=eig(S);  % the principal stresses are the eigen values of the tensor S;

% Method of Maxpscheck. 11/19/2019
% current stress tensor [sxx,sxy; sxy, syy]
S=[s(1), s(3); s(3), s(2)];
% rotation angle
gtheta=0.5*atan(2*s(3,:)./(s(1,:)-s(2,:)));    % in radians
% drop minor numerical error to avoid accumulating deflection 11/20/19
if abs(gtheta)<1e-9
    gtheta=0;
elseif abs(gtheta)>pi/2
    warning('The crack deflection can be overwhelming')
    gtheta=1/3*gtheta;
end
A=[cos(gtheta), sin(gtheta);-sin(gtheta), cos(gtheta)];
% transformed stress tensor
St=A*S*A';
% principal stress along x' axis
psx=St(1,1);
% principal stress along y' axis
psy=St(2,2);
psmax=max(psx,psy);
% if tip==1 % the right tip is relatively easy to manage
% %     if psx>=psy
% %         % psmax is psx rotated by gtheta from x axis
% %         % then the crack direction should be normal to x', then y
% % %         gtheta=gtheta+sign(omega)*pi/2;
% % 		warning('the maximum principal stress direction has changed')
% %         % else
% %         %     % psmax is psy rotated by gtheta from y axis
% %         %     % then the crack direction should be normal to y'
% %     end
% else % if tip ==2, the left tip
% %     if psx>=psy
% %         % psmax is psx rotated by gtheta from x axis
% %         % then the crack direction should be normal to x', then y
% % %         gtheta=gtheta+sign(omega)*pi/2;
% % 		warning('the maximum principal stress direction has changed')
% %         % else
% %         %     % psmax is psy rotated by gtheta from y axis
% %         %     % then the crack direction should be normal to y'
% %     else
if tip==2
        gtheta=gtheta+pi;
end

obj.PVariable=[obj.PVariable,psmax];
obj.Theta=[obj.Theta,gtheta];
end

