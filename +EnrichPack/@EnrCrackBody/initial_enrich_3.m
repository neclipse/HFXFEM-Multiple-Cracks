function initial_enrich_3(obj,varargin)
% method of EnrCrackBody
% enrich elems using all enrichment functions. 
for iE=1:length(obj.NewElems)
    for ienf=1:length(obj.Myenfs)
        myenf=obj.Myenfs{ienf};
        myenf.enrichelem(obj.Elemdict(obj.NewElems(iE)),obj.Id);
    end
end
end

