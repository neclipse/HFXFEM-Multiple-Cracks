function update_enrich_1(obj,varargin)
% method of EnrCrackBody update geometry and associated elements
% and nodes, change the enrich flag of these elements and nodes, but do not
% update the subdomain and update the enrichment function.
%% Really grow the crack if cutflag is not true for any tip
nodes=NaN(1,100);
elems=NaN(1,100);
inode=1;ielem=1;
for itip=1:length(obj.Mytips)
    % really propagate mygeo if needed
    obj.Mytips(itip).realgrow;
    nodes(inode:inode+length(obj.Mytips(itip).NewNodes)-1)=obj.Mytips(itip).NewNodes;
    elems(ielem:ielem+length(obj.Mytips(itip).NewElems)-1)=obj.Mytips(itip).NewElems;
    inode=inode+length(obj.Mytips(itip).NewNodes);
    ielem=ielem+length(obj.Mytips(itip).NewElems);
end
nodes=nodes(1:inode-1);
elems=elems(1:ielem-1);
obj.NewNodes=nodes;
obj.NewElems=elems;
%% Really update crackbody object info and the enrichment
if any([obj.Mytips.Growcheck.Growflag])
    display(growflags,'crack growth');
    perforated=false;
    %% update basic info
    obj.setinteractedelem;
    obj.setenrichednode;
    for ienf=1:length(obj.Myenfs)
        obj.Myenfs{ienf}.Lsv=obj.Mygeo.Phi;
    end
    %% update the enrichement
    % nodes then elems, unnecessary for now as std nodes are also
    % enriched for additional dofs, only they are kept zeros
    for iN=1:length(nodes)
        obj.Nodedict(nodes(iN)).setenrich(obj.Id);
    end
    for iE=1:length(elems)
        % set enrich flag and find the standard nodes within the
        % element
        obj.Elemdict(elems(iE)).setenrich(obj.Id);
        % divide the element into triangular subdomains for 2d
        % integral, this module is moved to update_enrich_2. 11/06/20
%         obj.Elemdict(elems(iE)).subdomain(obj.Id);
        % find the gaussian points on the crack for line integral,
        % p=3 to make the gauss quadrature accurate enough for the
        % line integral. (not sure if it is really useful, 03122019)
        % Change p=2 (defautl value) on 06072019.
        obj.Elemdict(elems(iE)).linegauss(obj.Id,obj.Cohesive,perforated);
    end
    %% update enrich all new nodes inside the enriched elements
    % the standard nodes have the enriched dofs but keep zero
    % no need to have these after 10/02/20 as udof==2, pdof==1,
    % dofs==3, always.
%     for iN=1:length(nodes)
%         node=nodes(iN);
%         for ienf=1:length(obj.Myenfs)
%             myenf=obj.Myenfs{ienf};
%             myenf.addnodedofs(obj.Nodedict(node),obj.Id);
%         end
%     end
    %% update enrich new elemdict,this module is update_enrich_2. 11/06/20
%     for iE=1:length(elems)
%         elem=elems(iE);
%         % use geometeric info to divide the elem into subdomains
%         % for integral purpose.
%         %mygeo.subdomain(elem);
%         for ienf=1:length(obj.Myenfs)
%             myenf=obj.Myenfs{ienf};
%             myenf.enrichelem(obj.Elemdict(elem),obj.Id);
%         end
%     end
end
obj.checkactive;
end

