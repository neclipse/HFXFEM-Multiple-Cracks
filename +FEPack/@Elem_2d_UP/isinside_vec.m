function [flagie,flagi,flage,flagoe,area] = isinside_vec(obj,plist)                    
% Check if point (x,y) lies inside this element or on some edge of the element
% can be replaced by the MATLAB built-in "inpolygon" function 03282019
p0=plist;														% a list of point(x,y)
n=size(p0,1);
area=zeros(n,4);									% array of area of four triangles
if isempty(obj.Area)
    obj.calarea;
end
if isempty(obj.Length)
    obj.callength;
end
p1=repmat([obj.X(1),obj.Y(1)],n,1);								% For vectorization
p2=repmat([obj.X(2),obj.Y(2)],n,1);
p3=repmat([obj.X(3),obj.Y(3)],n,1);
p4=repmat([obj.X(4),obj.Y(4)],n,1);
l1=repmat(obj.Length(1),n,1);									% element length 1
l2=repmat(obj.Length(2),n,1);
l3=repmat(obj.Length(3),n,1);
l4=repmat(obj.Length(4),n,1);
dist1=real(sqrt((p0(:,1)-p1(:,1)).^2+(p0(:,2)-p1(:,2)).^2));												% distance from p0 to p1;
dist2=real(sqrt((p0(:,1)-p2(:,1)).^2+(p0(:,2)-p2(:,2)).^2));
dist3=real(sqrt((p0(:,1)-p3(:,1)).^2+(p0(:,2)-p3(:,2)).^2));
dist4=real(sqrt((p0(:,1)-p4(:,1)).^2+(p0(:,2)-p4(:,2)).^2));
% Areas of 4 triangles
s1=(dist1+dist2+l1)./2;
area(:,1)=real(sqrt(s1.*(s1-dist1).*(s1-dist2).*(s1-l1)));
s2=(dist2+dist3+l2)./2;
area(:,2)=real(sqrt(s2.*(s2-dist2).*(s2-dist3).*(s2-l2)));
s3=(dist3+dist4+l3)./2;
area(:,3)=real(sqrt(s3.*(s3-dist3).*(s3-dist4).*(s3-l3)));
s4=(dist4+dist1+l4)./2;
area(:,4)=real(sqrt(s4.*(s4-dist4).*(s4-dist1).*(s4-l4)));
% Compare the sum area to the Element area to see if the point is inside or on the edge of the element
sumarea=sum(area,2);										% sum all four areas along the row
err=(sumarea-obj.Area)/obj.Area;
flagie=(err<1e-7);                                              % True means inside or on the edge, False means outside
% If not all areas on the same row is nonzero, then the point is on the edge of the element
% flagnode=(area<0.005*obj.Area);
% flagnode=sum(flagnode,2);
% flagnode=flagnode==2;                                           % if two areas is almost zero, then it means the point is too close to a node
flage=~all(area>1e-8,2);										% True means on the edge line, False means not
flagi=flagie & ~flage;                                          % true means inside the element and not on the edges and not too close to a node
flagoe=~flagie & flage;											% true means outside the element but on the edge line
flage(flagoe)=false;                                            % reset to false when the point is actually on the extension of the edge
end
