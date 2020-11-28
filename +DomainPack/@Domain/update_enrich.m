function update_enrich(obj)
%UPDATE_ENRICH of Domain class, moved back from Newton_Raphson class
%11/27/20 as this update comes after the coverged solution.
%Higher level function to guide all enrichitems

if ~isempty(obj.EnrichItems)
    %1. do the postprocessing for the current enrichitems before growing the
%new increment.
    for ienr=1:length(obj.EnrichItems)
        enritem=obj.EnrichItems{ienr};
        enritem.NewElems=[];
        enritem.NewNodes=[];
        enritem.postprocess(obj.NewtonRaphson.Dt);      % obtain crack aperture and other practical information.
    end
    
    %2. Really grow the cracks and update enrich for every enritem on level 1:
    %found the newelems and set enrich, apply the linegauss
    for ienr=1:length(obj.EnrichItems)
        obj.EnrichItems{ienr}.update_enrich_1;
    end
    % 3. Update the enritem on level 2: apply subdomain and update
    % enf.enrichelem after all enichitems have finished level 1
    for ienr=1:length(obj.EnrichItems)
        obj.EnrichItems{ienr}.update_enrich_2;
    end
end
    % The updatedofarray_enriched should be run after the obj.storage.
    % 11/27/2020.

end


