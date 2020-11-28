function initiate_enrich(obj)
%initiate_ENRICH function of Domain class to initiate the enrichitems
% 11/05/2020
if ~isempty(obj.EnrichItems)
    for ienr=1:length(obj.EnrichItems)
        mygeo=obj.EnrichItems{ienr}.Mygeo;
        enfh=EnrichPack.EnFHeaviside(mygeo.Minelength,mygeo.Phi); % smoothed, type 1 heaviside is default
        enfrd=EnrichPack.EnFRidge(mygeo.Phi);
        % May add another junction enrichment function for minor crack tip
        obj.EnrichItems{ienr}.Myenfs={enfh,enfrd};
        obj.EnrichItems{ienr}.initial_enrich_1;
    end
    % 2. Level 2: generate the line gaussian points for all interacted
    % elems after level 1: the interacted elems are properly set enriched
    % with correct EnrichNum, also generate the domain gaussian points.
    for ienr=1:length(obj.EnrichItems)
        obj.EnrichItems{ienr}.initial_enrich_2;
    end
    % 3. enf.enrichelem after all enichitems have finished level 2
    for ienr=1:length(obj.EnrichItems)
        obj.EnrichItems{ienr}.initial_enrich_3;
    end
    % update the whole dof array
    obj.updatedofarray_enriched;  % update the linsystem inside
end
end

