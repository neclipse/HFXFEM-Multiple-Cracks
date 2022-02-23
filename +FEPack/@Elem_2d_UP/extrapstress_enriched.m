function [p,stressp]=extrapstress_enriched(obj,xi,eta)
%extrapstress_enriched: function of elem2d, will be called in
%enrichitem.postprocess to estimate the stress at other line gausian points
xg=[obj.EnrichGauss.Xi];
yg=[obj.EnrichGauss.Eta];
dist=sqrt((xg-xi).^2+(yg-eta).^2);
[~,I]=sort(dist);   % assended order
p=(obj.EnrichGauss(I(1)).P+...
    obj.EnrichGauss(I(2)).P+...
    obj.EnrichGauss(I(3)).P)/3;
stressp=(obj.EnrichGauss(I(1)).Stressp+...
    obj.EnrichGauss(I(2)).Stressp+...
    obj.EnrichGauss(I(3)).Stressp)/3;
% column stress vector: [rep;sx;sy;sxy;sz]
stressp=[stressp(1),stressp(3);stressp(3),stressp(2)]; % convert stressp to tensor form
end


