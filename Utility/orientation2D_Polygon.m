function [ c ] = orientation2D_Polygon( V )
%This algorithm is to test the orientation of an arbitrary cyclic polygon
% c is positive when the polygon is oriented counterclockwise, vice versa.
%   Find the rightmost lowest vertex and then calculate cross product of
%   edges fore and aft of it.
n=size(V,1);
V(n+1,:)=V(1,:);
kmin=1;
xmin=V(kmin,1);
ymin=V(kmin,2);
%% Find the rightmost lowest node
for i=2:n
    if V(i,2)<=ymin
        if V(i,1)>xmin
            kmin=i;
            xmin=V(i,1);
            ymin=V(i,2);
        end
    end 
end
%% Calculate the cross product of two vectors at the node
if kmin==1;
    c=(V(1,1)-V(n,1))*(V(2,2)-V(n,2))-(V(2,1)-V(n,1))*(V(1,2)-V(n,2));
else
    c=(V(kmin,1)-V(kmin-1,1))*(V(kmin+1,2)-V(kmin-1,2))-(V(kmin+1,1)-V(kmin-1,1))*(V(kmin,2)-V(kmin-1,2));
end

end

