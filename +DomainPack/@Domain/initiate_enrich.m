function initiate_enrich(obj)
%initiate_ENRICH function of Domain class to initiate the enrichitems
% 11/05/2020
for ienr=1:obj.length(obj.EnrichItems)
    mygeo=obj.EnrichItems{ienr}.Mygeo;
    enfh=EnrichPack.EnFHeaviside(mygeo.Minelength,mygeo.Phi); % smoothed, type 1 heaviside is default
    enfrd=EnrichPack.EnFRidge(mygeo.Phi);
    % May add another junction enrichment function for minor crack tip
    obj.EnrichItems{ienr}.Myenfs={enfh,enfrd};
    obj.EnrichItems{ienr}.initial_enrich;
end
% update the whole dof array
obj.updatedofarray_enriched;  % update the linsystem inside
end

