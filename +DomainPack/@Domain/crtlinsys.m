function obj=crtlinsys(obj)
%CRTLINSYS Summary of this function goes here
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
linsys=DomainPack.LinSysCreater(totdofs,totudofs,totpdofs,obj.ElemType,obj.MatType,...
    obj.ElemDict,obj.NodeDict,obj.Preprocess.BCTable_line,obj.Preprocess.BCTable_imbalance,...
    obj.Preprocess.BCTable_stress,obj.Preprocess.BCTable_traction,...
    obj.Preprocess.BCTable_node,obj.Preprocess.BCTable_disp,...
    obj.Preprocess.BCTable_flow,obj.EnrichItems);
obj.LinSysCrt=linsys;
% Added on 04232019
obj.assign_arbitrary_flow;
end

