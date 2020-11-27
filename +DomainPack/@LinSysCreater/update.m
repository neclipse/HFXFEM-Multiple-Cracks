function  update(obj,totdofs,totudofs,totpdofs,elemdict,nodedict,enrichitems)
%UPDATE The core of Domain.updatelinsys
ctotudofs=totudofs-obj.Nudof;
ctotpdofs=totpdofs-obj.Npdof;
%   Update the linsys due to the change of dimension
obj.Drow=totdofs;               % All mixed dofs
obj.Dcol=totdofs;
obj.Nudof=totudofs;             % All displacement dofs including the enriched
obj.Npdof=totpdofs;             % All pressure dofs including the enriched
% elemdict and nodedict should not be changed if the mesh is not adaptively
% updated 
obj.ElemDict=elemdict;       % Handles of all elements
obj.NodeDict=nodedict;       % Handles of all nodes
% if no new enrichitem is created
obj.EnrichItems=enrichitems;
obj.UO=[obj.UO;zeros(ctotudofs,1)];
obj.UOt1=[obj.UOt1;zeros(ctotudofs,1)];
obj.UOt2=[obj.UOt2;zeros(ctotudofs,1)];
obj.PO=[obj.PO;zeros(ctotpdofs,1)];
obj.POt1=[obj.POt1;zeros(ctotpdofs,1)];
end

