function Rtipelems=Toedge2( obj, varargin)
%TOEDGE2 To account for the tips initially on the element edge
%Called in obj.initiate or obj.update to achieve the following functions:
% Find these nodes on the tip edges
% 2(update on 02262019): tag the detected nodes as standard nodes
if ~isempty(varargin)
    Rtips=varargin{1};
else
    Rtips=obj.Rtips;
end
% Loop over every real tips
Rtipelems=zeros(length(Rtips),1);
for itip=1:length(Rtips)
    tip_ini=obj.Tips(Rtips(itip),:);
    rtipelems=obj.findelems(tip_ini,'in_edge');  % different from obj.initiate
    if length(rtipelems) ==2||4            % tip initially on an edges or an node
        % however, obj.Rtipelements=[]; the possible element containing the
        % tip is stored in obj.Intelements
        % we need find the element containing the tip
        rtipelem=setdiff(rtipelems,obj.Intelements);
        rtipelem=setdiff(rtipelems,rtipelem);
        Rtipelems(itip)=rtipelem;
        elem=obj.Elemdict(rtipelem);
        % find on which edge the tip lies
        [~,~,~,~,area] = isinside_vec(elem,tip_ini);
        [~,edge_ind]=min(area);
        % remove the two nodes on the found edge from the obj.Nodes
        if edge_ind<length(elem.NodList)
            to_remove=elem.NodList([edge_ind,edge_ind+1]);
        elseif edge_ind==length(elem.NodList)
            to_remove=elem.NodList([edge_ind,1]);
        end
        obj.Stdnodes(2*Rtips(itip)-1:2*Rtips(itip))=to_remove;
%         obj.Nodes=setdiff(obj.Nodes,to_remove);
%         for j=1:length(to_remove)
%             obj.Phi(obj.Phi(:,1)==to_remove(j),:)=[];
%         end
    end
end
%     obj.Stdnodes=unique(obj.Stdnodes);
end
