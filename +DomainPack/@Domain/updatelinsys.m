function obj=updatelinsys(obj)
%update linsys here as the 
%   Detailed explanation goes here
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
obj.LinSysCrt.update(totdofs,totudofs,totpdofs,obj.ElemDict,obj.NodeDict,obj.EnrichItems)

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

