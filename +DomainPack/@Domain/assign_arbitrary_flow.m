function obj=assign_arbitrary_flow( obj )
%ASSIGN_ARBITRARY_FLOW Summary of this function goes here (04232019)
% Method of the Domain class
% After the elements have been substantiated, update the
% Linsys.BCTable_flow using Preprocess.QTable_node if it is not empty.
%   Detailed explanation goes here
if ~isempty(obj.Preprocess.QTable_node)
    QTable=obj.Preprocess.QTable_node;
    BCTable_flow=zeros(1000,3);
    iedge=1;
    for iline=1:size(QTable,1)
        localnodes=QTable{iline,1}; % local index of the nodes along the edge
        points=QTable{iline,2};
        q=QTable{iline,3};
        mesh=obj.Preprocess.Mesher;
        elems=mesh.findelems(points,obj.ElemDict,'inside');
        % Find the edge using the localnodes
        for ielem=1:length(elems)
            nodes=obj.ElemDict(elems(ielem)).NodList(localnodes);
            BCTable_flow(iedge,:)=[nodes,q];
            iedge=iedge+1;
        end
    end
    BCTable_flow=BCTable_flow(1:iedge-1,:);
    % Append the new flow table to the linsys.bctableq, which will be used
    % in the linsys.crtRHS_UP automatically.
    obj.LinSysCrt.BCTableq=[obj.LinSysCrt.BCTableq;{BCTable_flow}];
end
end

