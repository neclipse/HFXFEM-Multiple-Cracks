function calelemstress(obj)
% Avearge the integration results to nodal values, NOT ACCURATE, ESPICICALLY ON THE BOUNDARIES
% collect the temporary stress results from elem.calstress and
% elem.calstress_enriched, the strains are not collected 05282019
% Need to ensure that the stress values stored in the nodes are cleared
% because the calculation need a clean slate. 07132019
for in=1:length(obj.NodeDict)
    obj.NodeDict(in).Stress=zeros(5,1) ;
    obj.NodeDict(in).Stressp=zeros(5,1) ;
end
% This is a wrong use of cell array "newelems=[obj.EnrichItems.NewElems];"
% obj.EnrichItems is a cell array not an object array
% newelems=[];
% for ienr=1:length(obj.EnrichItems)
%     newelems=[newelems,obj.EnrichItems{ienr}.NewElems];
% end
% At a second thought, change cell array to object array.
newelems=unique([obj.EnrichItems.NewElems]);
% use storage flag to indicate that nodal stress is to be accumulated.
% 12/22/20
storage=true;
for ie=1:length(obj.ElemDict)
    elem=obj.ElemDict(ie);
    % Note that the newly enriched element does not have valid enriched dofs
    % Cannot use calstress_enriched to update stress. 09272019
%     try
        if any(elem.RealEnrich) && ~any(newelems==elem.Ind)     % Enriched elements
            %         ienrich= find(elem.Enrich);
            elem.calstress_enriched(storage);
        else                    % Standard elements
            elem.calstress(storage);
        end
%     catch
%         fprintf('Problem to call calstress_enriched for No.% element',elem.Ind);
%         continue;
%     end
end
end