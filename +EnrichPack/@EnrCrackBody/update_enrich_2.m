function update_enrich_2(obj,varargin)
% method of EnrCrackBody
% update the subdomain and line gaussian points
growchecks=[obj.Mytips.Growcheck];
if any([growchecks.Growflag])
    perforated=false;
    for iE=1:length(obj.NewElems)
        % divide the element into triangular subdomains for 2d integral
        % THINK IF THE SUBDOMAIN AND ERNICHELEME CAN BE CALLED HERE. 11/06/20
        % Yes, subdomain can be called here, together with linegauss.
        % 11/27/20
        obj.Elemdict(obj.NewElems(iE)).subdomain;
        obj.Elemdict(obj.NewElems(iE)).linegauss(obj.Id,obj.Cohesive,perforated,obj.Alpha);
    end
end
end

