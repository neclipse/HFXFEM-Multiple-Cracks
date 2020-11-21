function initial_enrich_2(obj,varargin)
% method of EnrCrackBody
% update the subdomain and update the enrichment function. 
for iE=1:length(obj.NewElems)
    % divide the element into triangular subdomains for 2d integral
    % THINK IF THE SUBDOMAIN AND ERNICHELEME CAN BE CALLED HERE. 11/06/20
    obj.Elemdict(obj.NewElems(iE)).subdomain;
    for ienf=1:length(obj.Myenfs)
        myenf=obj.Myenfs{ienf};
        myenf.enrichelem(obj.Elemdict(obj.NewElems(iE)),obj.Id);
    end
end
end

