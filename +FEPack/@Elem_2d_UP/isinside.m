function flag = isinside(obj,x,y)                                   % Check if point (x,y) is inside this element
		p1 = [x,y];
		temp = 0;
		for i = 1 : obj.NoNodes                                         % calculate the area of the triangle made by the point and every side of the element
			p2= [obj.X(i),obj.Y(i)];
			if i< obj.NoNodes
				p3 = [obj.X(i+1), obj.Y(i+1)];
			else
				p3 = [obj.X(1),obj.Y(1)];
			end
			l1 = norm(p1-p2); l2 = norm(p2-p3); l3 = norm (p3-p1);
			s = (l1+l2+l3)/2;
			temp=temp+sqrt(s*(s-l1)*(s-l2)*(s-l3));
		end
		if isempty(obj.Area)
			obj.calarea;
		end
		if temp == obj.Area
			flag = true;                                                % point (x,y) is inside the element
		else
			flag = false;                                               % point(x,y) is outside the element
		end
end
