function crtpsdbound( obj )
%CRTPSDBOUND Summary of this function goes here
%   Detailed explanation goes here
nodes=cell(size(obj.Psdboundind,1),1);
for ipb=1:size(obj.Psdboundind,1)
    ir=logical(obj.BCTable_line(:,1)==obj.Psdboundind(ipb,1));
    edgetable=obj.BCTable_line(ir,2:3);                 % edges on this psdbound
    edgevec=edgetable(:);                               % convert matrix to vector
    nodes{ipb}=unique(edgevec,'sorted');                % all nodes on this psdbound
end
obj.Psdnodes=nodes;
end

