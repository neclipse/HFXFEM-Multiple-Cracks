function initiate_enrich(obj)
%initiate_ENRICH function of Domain class to initiate the enrichitems
% 11/05/2020
if ~isempty(obj.EnrichItems)
    for ienr=1:length(obj.EnrichItems)
        mygeo=obj.EnrichItems(ienr).Mygeo;
        enfh=EnrichPack.EnFHeaviside(mygeo.Minelength,mygeo.Phi); % smoothed, type 1 heaviside is default
        enfrd=EnrichPack.EnFRidge(mygeo.Phi);
        % May add another junction enrichment function for minor crack tip
        obj.EnrichItems(ienr).Myenfs={enfh,enfrd};
        obj.EnrichItems(ienr).initial_enrich_1;
    end
    % 2. Level 2: generate the line gaussian points for all interacted
    % elems after level 1: the interacted elems are properly set enriched
    % with correct EnrichNum, also generate the domain gaussian points.
    for ienr=1:length(obj.EnrichItems)
        obj.EnrichItems(ienr).initial_enrich_2;
    end
    % 3. enf.enrichelem after all enichitems have finished level 2
%     for ienr=1:length(obj.EnrichItems)
%         obj.EnrichItems(ienr).initial_enrich_3;
%     end

    % Alternative approach to fix issue #19. The old approach should have
    % no problems but it is changed just to be consistent with
    % Domain.update_enrich 12/25/2020.
    % NewElems already excludes the elements that only involve with smeared
    % cracks. 02/04/2021.
    allnewelems=unique([obj.EnrichItems.NewElems]);
    % Loop from the newelem level
    for ielem=1:length(allnewelems)
        newelem=obj.ElemDict(allnewelems(ielem));
        % loop over all real EnrichItems involved with the newelem
        for i=1:newelem.EnrichNum
            ienr=newelem.get_realenrichind(i); % 
            EnrItem=obj.EnrichItems(newelem.Enrich(ienr));
            % loop over the Myenfs of the EnrItem
            for ienf=1:length(EnrItem.Myenfs)
                myenf=EnrItem.Myenfs{ienf};
                myenf.enrichelem(newelem,EnrItem.Id);
            end
        end
    end
    % move updatedof_enriched outside to domain.running directly to avoid
    % the return of obj. 12/04/2020.
end
end

