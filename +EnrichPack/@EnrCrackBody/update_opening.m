function update_opening(obj) 
% concrete method of EnrCrackBody

% To get the Heaviside type
for ienf=1:length(obj.Myenfs)
    enf=obj.Myenfs{ienf};
    if isa(enf,'EnrichPack.EnFHeaviside')
        enfh=enf;
        break;
    end
end
type=enfh.Type;
for iE=1:length(obj.INTELEM)
    elem=obj.INTELEM(iE);
    Jacob=elem.JacobianMatDict(obj.Id);
	% Enriched Dofs
	ue=Jacob.Un2iEnr;
    linegauss=elem.LineGaussDict{obj.Id};
   for ig=1:length(lg)
       lg=linegauss(ig);
       % update the crack opening using enriched udofs
       lg=lg.updatecrackopening(ue,type);
   end
end
end
