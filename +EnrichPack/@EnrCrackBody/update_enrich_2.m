function update_enrich_2(obj,varargin)
% method of EnrCrackBody
% update the subdomain and line gaussian points for obj.TransElems
for iE=1:length(obj.TransElems)
    % divide the element into triangular subdomains for 2d integral
    % THINK IF THE SUBDOMAIN AND ERNICHELEME CAN BE CALLED HERE. 11/06/20
    % Yes, subdomain can be called here, together with linegauss.
    % 11/27/20
    obj.Elemdict(obj.TransElems(iE)).subdomain;
    ind=obj.Elemdict(obj.TransElems(iE)).Enrich==obj.Id;
    % LineGauss already called for previously smeared elements.
    linegaussdict=obj.Elemdict(obj.TransElems(iE)).LineGaussDict{ind};
    for ip=1:length(linegaussdict)
        linegaussdict(ip)=linegaussdict(ip).transit;
    end
    obj.Elemdict(obj.TransElems(iE)).LineGaussDict{ind}=linegaussdict;
end
% create the subdomain and line gaussian points for the rest of obj.NewElems
newpropagatedelems=setdiff(obj.NewElems,obj.TransElems);
initialmode=5; % newly propagated elems starts with initial tensile mode
for iE=1:length(newpropagatedelems)
    % divide the element into triangular subdomains for 2d integral
    % THINK IF THE SUBDOMAIN AND ERNICHELEME CAN BE CALLED HERE. 11/06/20
    % Yes, subdomain can be called here, together with linegauss.
    % 11/27/20
    obj.Elemdict(newpropagatedelems(iE)).subdomain;
    % But: for newly created crack, should not use obj.Alpha but the
    % default value. Issue #20.
    obj.Elemdict(newpropagatedelems(iE)).linegauss(obj.Id,obj.Cohesive,initialmode);
end
end

