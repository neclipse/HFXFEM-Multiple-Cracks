function initial_enrich_2(obj,varargin)
% method of EnrCrackBody
% update the subdomain and line gaussian points
for iE=1:length(obj.Interactedelem)
    % divide the element into triangular subdomains for 2d integral
    % THINK IF THE SUBDOMAIN AND ERNICHELEME CAN BE CALLED HERE. 11/06/20
    % Yes, subdomain can be called here, together with linegauss.
    % 11/27/20
    if obj.InitialMode~=2 
        % No need to run subdomain for smeared crack, if there is smeared
        % carck and open crack, the open crack will have the subdomain
        % called with no problem. 02/04/2021.
        obj.Elemdict(obj.Interactedelem(iE)).subdomain;
    end
    obj.Elemdict(obj.Interactedelem(iE)).linegauss(obj.Id,obj.Cohesive,obj.InitialMode,obj.Alpha);
end
end

