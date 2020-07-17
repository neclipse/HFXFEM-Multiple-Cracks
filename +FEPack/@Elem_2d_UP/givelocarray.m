function  givelocarray(obj,varargin)
%GIVELOCARRAY Calculate the global index of DOf
% if nargin==1 % all nodes are involved
%     nnodes=obj.NoNodes;
% elseif nargin==2 % only selected nodes are involved
%     nnodes=varargin{1};
% end
locarray=[obj.NodDict(:).DofArray];
locarrayU=[obj.NodDict(:).UArray];
locarrayP=[obj.NodDict(:).PArray];
obj.Locarray=locarray;
obj.LocarrayU=locarrayU;
obj.LocarrayP=locarrayP;
end



