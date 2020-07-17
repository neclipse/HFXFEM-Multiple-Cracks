function  givelocarray_enriched(obj,crackid,varargin)
%GIVELOCARRAY Calculate the global index of DOf
% if nargin==1 % all nodes are involved
%     nnodes=obj.NoNodes;
% elseif nargin==2 % only selected nodes are involved
%     nnodes=varargin{1};
% end
locarrayenr=[obj.NodDict(1).EnrDofArray{crackid},obj.NodDict(2).EnrDofArray{crackid},...
    obj.NodDict(3).EnrDofArray{crackid},obj.NodDict(4).EnrDofArray{crackid}];
locarrayUenr=[obj.NodDict(1).UenrArray{crackid},obj.NodDict(2).UenrArray{crackid},...
    obj.NodDict(3).UenrArray{crackid},obj.NodDict(4).UenrArray{crackid}];
locarrayPenr=[obj.NodDict(1).PenrArray{crackid},obj.NodDict(2).PenrArray{crackid},...
    obj.NodDict(3).PenrArray{crackid},obj.NodDict(4).PenrArray{crackid}];
% Directly store the locarray to the specified crack JacobianMat
obj.JacobianMatDict(crackid).LocarrayEnr=locarrayenr;
obj.JacobianMatDict(crackid).LocarrayUEnr=locarrayUenr;
obj.JacobianMatDict(crackid).LocarrayPEnr=locarrayPenr;
end



