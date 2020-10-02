function  givelocarray_enriched(obj,crackid,varargin)
%GIVELOCARRAY Calculate the global index of DOf
% if nargin==1 % all nodes are involved
%     nnodes=obj.NoNodes;
% elseif nargin==2 % only selected nodes are involved
%     nnodes=varargin{1};
% end
% Now we need to use logical indexing to retrieve the exact index k where the
% requested crackid is stored in each of the nodes.
k1=obj.NodDict(1).Enrich==crackid;
k2=obj.NodDict(2).Enrich==crackid;
k3=obj.NodDict(3).Enrich==crackid;
k4=obj.NodDict(4).Enrich==crackid;
locarrayenr=[obj.NodDict(1).EnrDofArray{k1},obj.NodDict(2).EnrDofArray{k2},...
    obj.NodDict(3).EnrDofArray{k3},obj.NodDict(4).EnrDofArray{k4}];
locarrayUenr=[obj.NodDict(1).UenrArray{k1},obj.NodDict(2).UenrArray{k2},...
    obj.NodDict(3).UenrArray{k3},obj.NodDict(4).UenrArray{k4}];
locarrayPenr=[obj.NodDict(1).PenrArray{k1},obj.NodDict(2).PenrArray{k2},...
    obj.NodDict(3).PenrArray{k3},obj.NodDict(4).PenrArray{k4}];
% Directly store the locarray to the specified crack JacobianMat
k=obj.Enrich==crackid;
obj.JacobianMatDict(k).LocarrayEnr=locarrayenr;
obj.JacobianMatDict(k).LocarrayUEnr=locarrayUenr;
obj.JacobianMatDict(k).LocarrayPEnr=locarrayPenr;
end



