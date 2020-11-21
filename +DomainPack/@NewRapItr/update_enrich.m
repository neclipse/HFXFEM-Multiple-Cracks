function update_enrich(obj)
%UPDATE_ENRICH of Newton_Raphson Iterator class
%Higher level function to guide all enrichitems
%1. do the postprocessing for the current enrichitems before growing the
%new increment.
if ~isempty(obj.LinSysCrt.EnrichItems)
    for ienr=1:length(obj.LinSysCrt.EnrichItems)
        enritem=obj.LinSysCrt.EnrichItems{ienr};
        enritem.NewElems=[];
        enritem.NewNodes=[];
        enritem.postprocess(obj.Dt);      % obtain crack aperture and other practical information.
    end
end
%2. Really grow the cracks and update enrich for every enritem on level 1:
%found the newelems and set enrich, apply the linegauss
for ienr=1:length(obj.LinSysCrt.EnrichItems)
    obj.LinSysCrt.EnrichItems{ienr}.update_enrich_1;
end
% 3. Update the enritem on level 2: apply subdomain and update
% enf.enrichelem after all enichitems have finished level 1
for ienr=1:length(obj.LinSysCrt.EnrichItems)
    obj.LinSysCrt.EnrichItems{ienr}.update_enrich_2;
end

