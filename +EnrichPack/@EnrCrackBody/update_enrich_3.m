function update_enrich_3(obj,varargin)
% method of EnrCrackBody
% use enrichment functions to enrich the elem properly.
growchecks=[obj.Mytips.Growcheck];
if any([growchecks.Growflag])
    for iE=1:length(obj.NewElems)
        for ienf=1:length(obj.Myenfs)
            myenf=obj.Myenfs{ienf};
            myenf.enrichelem(obj.Elemdict(obj.NewElems(iE)),obj.Id);
        end
    end
end
end

