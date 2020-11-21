function area = calarea (obj)  
% Calculate element area
% if mod(obj.NoNodes,3) == 0
%     % use side lengths to calculate triangle area, Heron's formula
%     if isempty(obj.Length)
%         obj.callength;
%     end
%     s= sum(obj.Lengths)/2;
%     area=sqrt(s*(s-obj.Length(1))*(s-obj.Length(2))*(s-obj.Length(3)));
% elseif mod(obj.NoNodes, 4) == 0
%     if isempty(obj.Length)
%         obj.callength;
%     end
%     diagonal=sqrt((obj.X(1)-obj.X(3))^2+(obj.Y(1)-obj.Y(3))^2);
%     s1=(obj.Length(1)+obj.Length(2)+diagonal)/2;
%     s2=(obj.Length(3)+obj.Length(4)+diagonal)/2;
%     area1 = sqrt(s1*(s1-obj.Length(1))*(s1-obj.Length(2))*(s1-diagonal));
%     area2 = sqrt(s2*(s2-obj.Length(3))*(s2-obj.Length(4))*(s2-diagonal));
%     area = area1+area2;
% end
% Use the builtin method polyarea 03292019
x=obj.X;
y=obj.Y;
area=polyarea(x,y); % polyarea can be used when x,y are right ordered (no intersection.)
obj.Area = area;
end
% function   calarea( obj )
% %CALAREA Method of Elem_2d_EP to calculate the area of the element
% x=obj.X;
% y=obj.Y;
% obj.Area=0.5*abs((x(2)-x(4))*(y(3)-y(1))-(x(3)-x(1))*(y(2)-y(4)));
% end

