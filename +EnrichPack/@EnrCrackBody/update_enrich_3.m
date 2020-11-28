function update_enrich_3(obj,varargin)
% method of EnrCrackBody
% use enrichment functions to enrich the elem properly.
if any([obj.Mytips.Growcheck.Growflag])   
    for iE=1:length(obj.NewElems)
        for ienf=1:length(obj.Myenfs)
            myenf=obj.Myenfs{ienf};
            myenf.enrichelem(obj.Elemdict(obj.NewElems(iE)),obj.Id);
        end
    end
end
end

