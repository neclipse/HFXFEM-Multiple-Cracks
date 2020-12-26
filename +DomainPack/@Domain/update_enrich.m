function update_enrich(obj)
%UPDATE_ENRICH of Domain class, moved back from Newton_Raphson class
%11/27/20 as this update comes after the coverged solution.
%Higher level function to guide all enrichitems

if ~isempty(obj.EnrichItems)
    %1. do the postprocessing for the current enrichitems before growing the
%new increment.
    for ienr=1:length(obj.EnrichItems)
        enritem=obj.EnrichItems(ienr);
        enritem.NewElems=[];
        enritem.NewNodes=[];
        enritem.postprocess(obj.NewtonRaphson.Dt);      % obtain crack aperture and other practical information.
    end
    
    %2. Really grow the cracks and update enrich for every enritem on level 1:
    %found the newelems and set enrich, apply the linegauss
    for ienr=1:length(obj.EnrichItems)
        obj.EnrichItems(ienr).update_enrich_1;
    end
    % 3. Level 2: generate the line gaussian points for all interacted
    % elems after level 1: the interacted elems are properly set enriched
    % with correct EnrichNum
    for ienr=1:length(obj.EnrichItems)
        obj.EnrichItems(ienr).update_enrich_2;
    end
%     % 4. Update the enrichitem on level 2:  update enf.enrichelem after all
%     % enichitems have finished level 1 and level 2
%     for ienr=1:length(obj.EnrichItems)
%         obj.EnrichItems(ienr).update_enrich_3;
%     end
    % Fixing a BUG in issue #19: the existing crack did not enrich the
    % element with a new crack crossing but the existing enrichitem did not
    % have the element as NewElem. 12/25/2020
    % New approach: fetch all newelems and start the loop from the new elems
    % but not from the enrichitems. Domain. initiate_enrich can follow the
    % old approach or this new approach.
    allnewelems=unique([obj.EnrichItems.NewElems]);
    % Loop from the newelem level
    for ielem=1:length(allnewelems)
        newelem=obj.ElemDict(allnewelems(ielem));
        % loop over all EnrichItems involved with the newelem
        for ienr=1:newelem.EnrichNum
            EnrItem=obj.EnrichItems(newelem.Enrich(ienr));
            % loop over the Myenfs of the EnrItem
            for ienf=1:length(EnrItem.Myenfs)
                myenf=EnrItem.Myenfs{ienf};
                myenf.enrichelem(newelem,EnrItem.Id);
            end
        end
    end
end
    % The updatedofarray_enriched should be run after the obj.storage.
    % 11/27/2020.

end


