function update_enrich_1(obj,varargin)
% method of EnrCrackBody update geometry and associated elements
% and nodes, change the enrich flag of these elements and nodes, but do not
% update the subdomain and update the enrichment function.
%% Really grow the crack if cutflag is not true for any tip
nodes=NaN(100,1);
elems=NaN(100,1);
inode=1;ielem=1;
for itip=1:length(obj.Mytips)
    tipelem=obj.Mytips(itip).INTELEM;
    % only call checkgrow when the tipelem is already open
    if ~tipelem.Smeared(tipelem.Enrich==obj.Id)
        % really propagate mygeo if needed
        obj.Mytips(itip).realgrow;
        nodes(inode:inode+length(obj.Mytips(itip).NewNodes)-1)=obj.Mytips(itip).NewNodes;
        elems(ielem:ielem+length(obj.Mytips(itip).NewElems)-1)=obj.Mytips(itip).NewElems;
        inode=inode+length(obj.Mytips(itip).NewNodes);
        ielem=ielem+length(obj.Mytips(itip).NewElems);
    end
end
nodes=nodes(1:inode-1);
elems=elems(1:ielem-1);
obj.NewNodes=[obj.NewNodes;nodes];
obj.NewElems=[obj.TransElems;elems];
%% Really update crackbody object info and the enrichment
% Matlab does not allow second-level vectorized property retrieving.
% Then we use [obj.Mytips.Growcheck] to obtain a intermediate growcheck
% array. And we do another one on it to retrieve growflag.
growchecks=[obj.Mytips.Growcheck];
if any(obj.NewElems)
    if any([growchecks.Growflag])
        fprintf('\n The %d has crack growth at tips',obj.Id);
        disp(growchecks.Growflag);
    else
        fprintf('Part of the smeared crack %d opens.\n',obj.Id);
    end
    %% update basic info
    obj.setinteractedelem;
    obj.setenrichednode;
    % update the Lsv for enfs here
    for ienf=1:length(obj.Myenfs)
        obj.Myenfs{ienf}.Lsv=obj.Mygeo.Phi;
    end
    %% update the enrichement
    % nodes then elems, unnecessary for now as std nodes are also
    % enriched for additional dofs, only they are kept zeros
    for iN=1:length(obj.NewNodes) % setenrich all obj.NewNodes
        obj.Nodedict(obj.NewNodes(iN)).setenrich(obj.Id);
    end
    % only setenrich tip propagated elems, the opened smeared elems are
    % already updated in postprocess.
    for iE=1:length(elems) 
        % set enrich flag 
        elem=obj.Elemdict(elems(iE));
        elem.setenrich(obj.Id,4); % newly propagated segement initialmode=4.
        % Should use obj.Mygeo.Phi because the Phi values are updated based
        % on the new crack geometry. 01/22/2021 
        % [~,pnts,localpnts] = obj.Mygeo.intersection(elem, obj.Mygeo.Nodespool, obj.Mygeo.Phipool(:,2));
        % Now the Phipool are also correctly updated in tip.realgrow, so
        % this call does not need to be changed.
        [~,pnts,localpnts] = obj.Mygeo.intersection(elem);
        id=elem.Enrich==obj.Id;
        elem.LocalInt{id}=localpnts;
        elem.GlobalInt{id}=pnts;
        % divide the element into triangular subdomains for 2d
        % integral, this module is moved to update_enrich_2. 11/06/20
%         obj.Elemdict(elems(iE)).subdomain(obj.Id);
        % find the gaussian points on the crack for line integral,
        % p=3 to make the gauss quadrature accurate enough for the
        % line integral. (not sure if it is really useful, 03122019)
        % Change p=2 (defautl value) on 06072019.
%         elem.linegauss(obj.Id,obj.Cohesive,perforated); also moved to
%         update_enrich_2. 11/27/2020
    end
end
obj.checkactive;
end

