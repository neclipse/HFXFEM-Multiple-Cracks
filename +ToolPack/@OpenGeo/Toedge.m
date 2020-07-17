function Rtipelems=Toedge( obj )
%TOEDGE Extend the tips to edge for initial crack or grown crack
%Called in obj.initiate or obj.update to achieve the following functions:
% 1. Extend the tips to the edge according to the specified direction
% 2. Delete the nodes to be enriched in obj.Rtipnodes and obj.Nodes
% 2(update on 02262019): tag the detected nodes as standard nodes

% Loop over every real tips
Rtipelems=zeros(length(obj.Rtips),1);
for itip=1:length(obj.Rtips)
    tip_ini=obj.Tips(obj.Rtips(itip),:);
    rtipelems=obj.findelems(tip_ini,'to_edge');  % different from obj.initiate
    % function to describe the end segment;
    ends=obj.Ends{obj.Rtips(itip)};               % 2*2 matrix stroing the coordinates of the (1) tip and the (2) last point
    tan_vec=ends(1,:)-ends(2,:);               % 
    k1=tan_vec(2)/tan_vec(1);
    if k1==Inf
       k1=1e9;                              % to avoid possible error  
    end
    f1=@(x) k1*(x-tip_ini(1))+tip_ini(2);
    if length(rtipelems) == 1               % tip initially inside one element
        Rtipelems(itip)=rtipelems;
        elem=obj.Elemdict(rtipelems);
        XA=elem.X;                          % x coordinates of the four nodes inside elem
        YA=elem.Y;
        XB=[XA(2:end);XA(1)];
        YB=[YA(2:end);YA(1)];
        ka=(YB-YA)./(XB-XA);                % the slope vector of all edges
        ka(ka==Inf)=1e9;
        ka(ka==-Inf)=-1e9;
        xs=(k1*ends(1,1)-ka.*XA+YA-ends(1,2))./(k1-ka);     % derived solution
        % note xs may be Inf, but the point will be filtered out using
        % isinside_vec
        ys=f1(xs);                          % the coordinates of possible intersections
        plist=[xs,ys];
    %   test if the possible intersections are inside the elem
        [~,~,flage,~,~] = isinside_vec(elem,plist);
        ind=find(flage);
        if length(ind)<2
           error('One edge of tip element may not be correctly located during Toedge') 
        end
    %   find the intersection ahead of the crack tip, not the behind
        for i=1:length(ind)
            pfind=plist(ind(i),:);
            intvec=pfind-ends(1,:);
            psgn=sign(tan_vec*intvec');     % take the sign of the product of two vectors
            if psgn ==1
                edge_ind=ind(i);            % the index of the desired edge ahead of crack tip
                % insert the found pfind as the new Rtip
                index=size(obj.Segments,1)+1;
                newseg=[index,pfind];
                if obj.Rtips(itip)==2
                    obj.Segments=[newseg;obj.Segments];
                elseif obj.Rtips(itip)==1
                    obj.Segments=[obj.Segments;newseg];
                end
                % remove the two nodes on the found edge from the obj.Nodes
                if edge_ind<length(elem.NodList)
                    to_remove=elem.NodList([edge_ind,edge_ind+1]);
                elseif edge_ind==length(elem.NodList)
                    to_remove=elem.NodList([edge_ind,1]);
                end
%                 obj.Nodes=setdiff(obj.Nodes,to_remove);
                obj.Stdnodes=[obj.Stdnodes;to_remove'];
            end
        end
    elseif length(rtipelems) ==2||4            % tip initially on an edges
        % however, obj.Rtipelements=[]; the possible element containing the
        % tip is stored in obj.Intelements
        % we need find the element containing the tip
        continue;
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
        obj.Stdnodes=[obj.Stdnodes;to_remove'];
%         obj.Nodes=setdiff(obj.Nodes,to_remove);
%         for j=1:length(to_remove)
%             obj.Phi(obj.Phi(:,1)==to_remove(j),:)=[];
%         end
    end
end
end
