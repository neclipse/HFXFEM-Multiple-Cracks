function  initialize( obj )
%INITIALIZE initialize dofsold AND use upconf to assign these values to
%elements
%   UO, UOt1, UOt2, PO, POt1
% Give all zeros for the initial value
global inipore inidisp;
obj.UO=ones(obj.Nudof,1)*inidisp;
obj.UOt1=zeros(obj.Nudof,1); 
obj.UOt2=zeros(obj.Nudof,1);
obj.PO=ones(obj.Npdof,1)*(inipore); 
obj.POt1=zeros(obj.Npdof,1);
% Assign initial total pore pressure values or total displacements if necessary 
if ~isempty(obj.BCTabled)
   for id=1:length(obj.BCTabled) 
       disptable=obj.BCTabled{id};
       nodelist=disptable(1,:);
       idof=disptable(2,1);
       if idof ==1
           locarray=[obj.NodeDict(nodelist).UArray];  % Retrieve the boundary condition, 1-x_disp; 2-y-disp, 3-p
           locarray=locarray(idof:2:end);
           obj.UO(locarray)=disptable(3,:);
       elseif idof ==2
           locarray=[obj.NodeDict(nodelist).UArray];  % Retrieve the boundary condition, 1-x_disp; 2-y-disp, 3-p
           locarray=locarray(idof:2:end);
           obj.UO(locarray)=disptable(3,:);
       elseif idof ==3
           locarray=[obj.NodeDict(nodelist).PArray];  % Retrieve the boundary condition, 1-x_disp; 2-y-disp, 3-p
           obj.PO(locarray)=disptable(3,:);
       end
   end
end
% calculate initial RHS using initialstress
obj.initialRHS;
end

