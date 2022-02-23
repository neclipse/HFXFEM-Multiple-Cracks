function initial_enrich_1(obj)                                         % initialize the enrichment for this EnrichItem
% assign basic enrichment info to the selected elements and
% nodes
elems=obj.Interactedelem; % need to find the intersections no matter what initialmode.
nodes=obj.NewNodes; % if initialmode==2, no newnodes, no need to enrich nodes yet.
%% set the enrich flag = true and enrich_id=Enrichitem.Id, first
% nodes then elems
for iN=1:length(nodes)
    obj.Nodedict(nodes(iN)).setenrich(obj.Id);
end
for iE=1:length(elems)
    % set enrich flag and find the standard nodes within the
    elem=obj.Elemdict(elems(iE));
    elem.setenrich(obj.Id,obj.InitialMode);
    [~,pnts,localpnts] = obj.Mygeo.intersection(elem);
    id=elem.Enrich==obj.Id;
    elem.LocalInt{id}=localpnts;
    elem.GlobalInt{id}=pnts;
    % divide the element into triangular subdomains for 2d integral, this
    % method should be called after all cracks have been initially
    % enriched, so does the linegauss and initial enrich elemdict
    % module below. 10/30/20
    %
    %               elem.subdomain;
    %               subdomain and the enf.enrichelem should be moved out as there
    %               can be multiple cracks in this element initially. A good way
    %               is do this after finishing all enrcrack.initialenrich 11/20/20
    
    % find the gaussian points on the crack for line integral,
    % p=3 to make the gauss quadrature accurate enough for the
    % line integral. (not sure if it is really useful, 03122019)
    % Change p=2 (defautl value) on 06072019.
%     elem.linegauss(obj.Id,obj.Cohesive,obj.Perforated,obj.Alpha);
end
end