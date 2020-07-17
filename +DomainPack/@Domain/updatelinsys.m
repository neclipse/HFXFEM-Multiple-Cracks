function obj=updatelinsys(obj)
%update linsys here as the 
%   Detailed explanation goes here
linsys=obj.LinSysCrt;
if isempty(obj.EnrichItems)
    totdofs=obj.NoDofs;
    totudofs=obj.NoUDofs;
    totpdofs=obj.NoPDofs;
else
    totdofs=obj.NoDofs+obj.NoEnrDofs;
    totudofs=obj.NoUDofs+obj.NoUenrDofs;
    totpdofs=obj.NoPDofs+obj.NoPenrDofs;
end
% Find the dofs change
% ctotdofs=totdofs-linsys.Drow;
ctotudofs=totudofs-linsys.Nudof;
ctotpdofs=totpdofs-linsys.Npdof;
% update the basics in linsys
linsys.Drow=totdofs;               % All mixed dofs
linsys.Dcol=totdofs;
linsys.Nudof=totudofs;             % All displacement dofs including the enriched
linsys.Npdof=totpdofs;             % All pressure dofs including the enriched
% elemdict and nodedict should not be changed if the mesh is not adaptively
% updated
linsys.ElemDict=obj.ElemDict;       % Handles of all elements
linsys.NodeDict=obj.NodeDict;       % Handles of all nodes
% if no new enrichitem is created
linsys.EnrichItems=obj.EnrichItems;
linsys.UO=[linsys.UO;zeros(ctotudofs,1)];
linsys.UOt1=[linsys.UOt1;zeros(ctotudofs,1)];
linsys.UOt2=[linsys.UOt2;zeros(ctotudofs,1)];
linsys.PO=[linsys.PO;zeros(ctotpdofs,1)];
linsys.POt1=[linsys.POt1;zeros(ctotpdofs,1)];

obj.LinSysCrt=linsys;
% % Added on 04232019 to handle arbitrary flow condition
% obj.assign_arbitrary_flow;
% linsys.BCTable_line=obj.Preprocess.BCTable_line;
% linsys.BCTableim=obj.Preprocess.BCTable_imbalance;
% linsys.BCTables=obj.Preprocess.BCTable_stress;
% linsys.BCTablet=obj.Preprocess.BCTable_traction;
% linsys.BCTablen=obj.Preprocess.BCTable_node;
% linsys.BCTabled=obj.Preprocess.BCTable_disp;
% linsys.BCTableq=obj.Preprocess.BCTable_flow;

end

