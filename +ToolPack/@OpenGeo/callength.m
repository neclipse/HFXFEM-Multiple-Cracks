function length = callength( obj )
%CALLENGTH call external utility function arclength to calculate the total
%length of the crack represented by the Segments data.
xv=obj.Segments(:,2);
yv=obj.Segments(:,3);
[length,~]=arclength(xv,yv);
obj.Length=length;  % crack length
end

