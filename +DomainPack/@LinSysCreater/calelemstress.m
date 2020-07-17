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
newelems=[obj.EnrichItems{:}.NewElems];
for ie=1:length(obj.ElemDict)
    elem=obj.ElemDict(ie);
    % Note that the newly enriched element does not have valid enriched dofs
    % Cannot use calstress_enriched to update stress. 09272019
    if any(elem.Enrich) && ~any(newelems==elem.Ind)     % Enriched elements
%         ienrich= find(elem.Enrich);
        elem.calstress_enriched;
    else                    % Standard elements
        elem.calstress;
    end
end
end