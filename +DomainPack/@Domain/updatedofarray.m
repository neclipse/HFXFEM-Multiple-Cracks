function  obj=updatedofarray( obj )
%UPDATEDOFS Method of Domain class to update both the standard and the
%enriched dofs of all nodes within the domain.
%   Standard dofs
totndof=0;
totudof=0;
totpdof=0;
nnodes=length(obj.NodeDict);
nelem=length(obj.ElemDict);
for inod=1:nnodes
    [totndof,totudof,totpdof]=obj.NodeDict(inod).givedofarray(totndof,totudof,totpdof);
end
obj.NoDofs=totndof;                                 % number of total standard dofs.
obj.NoUDofs=totudof;
obj.NoPDofs=totpdof;
for ielem=1:nelem
    obj.ElemDict(ielem).givelocarray;
end
% crtpsddof using preprocess class info
if isempty(obj.PsdDofs)
    obj.Preprocess.crtpsdbound;
    psdnodes=obj.Preprocess.Psdnodes;
    psdind=obj.Preprocess.Psdboundind;
    psddofs=[];
    for ipb=1:size(psdind)
        nodesipb=psdnodes{ipb};
        idof=psdind(ipb,2);
        alldofs=[obj.NodeDict(nodesipb).DofArray];      % standard dofs of all found nodes on the boundary
        dofs=alldofs(idof:3:end);
        psddofs=[psddofs,dofs];
    end

    %% add the sparse nodes from preprocess.bctablen (03012019)
    bctablen=obj.Preprocess.BCTable_node;
    if ~isempty(bctablen)
        for inode=1:size(bctablen,1)
            idof=bctablen(inode,1);
            nodeid=bctablen(inode,2);
            dof=obj.NodeDict(nodeid).DofArray(idof);
            psddofs=[psddofs,dof];
        end
    end
    obj.PsdDofs=unique(psddofs);    %03/27/2019
end
% obj=obj.crtlinsys;                                      % create the linsystem here (03032019)
end