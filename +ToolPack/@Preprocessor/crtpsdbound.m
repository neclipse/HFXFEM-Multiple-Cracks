function crtpsdbound( obj )
%CRTPSDBOUND % Dirichlet boundary, in this case all bounaries are fixed
%   find the nodes on all Dirichlet boundaries, then the converger will
% not take the residual load vector on these nodes as true error but the 
% reactional force vector. Note the reactial force vector are only reported
% in the Domain.storage class for Postprocessor. 091120
nodes=cell(size(obj.Psdboundind,1),1);
for ipb=1:size(obj.Psdboundind,1)
    ir=logical(obj.BCTable_line(:,1)==obj.Psdboundind(ipb,1));
    edgetable=obj.BCTable_line(ir,2:3);                 % edges on this psdbound
    edgevec=edgetable(:);                               % convert matrix to vector
    nodes{ipb}=unique(edgevec,'sorted');                % all nodes on this psdbound
end
obj.Psdnodes=nodes;
end

