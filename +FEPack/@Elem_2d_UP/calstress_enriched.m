function calstress_enriched( obj)
%Calstress_enriched: function of elem2d, will be called in domain.storage
%Extrapolate/Average stress from subtriangles gaussian points to the nodes
%1. Loop over each node to find the closest two gaussian points
%2. Obtain the average stresses as the contribution of stresses from 
% this element to the node.
% obtain the coordinates of the enriched gaussian points
xg=[obj.EnrichGaussDict{1}.X];
yg=[obj.EnrichGaussDict{1}.Y];
% if isempty(obj.EnrichGaussDict{1}(1).Stressp)
    us=obj.Un2i;      % elemental iterative total displacement vector, u_sup(n+1)_sub(i+1)
    ps=obj.Pn2i;      % elemental iterative total pore pressure vector, p_sup(n+1)_sub(i+1)
    Jacob=obj.JacobianMatDict(1);
    ue=Jacob.Un2iEnr;
	pe=Jacob.Pn2iEnr;  
    for ig=1:length(obj.EnrichGaussDict{1})
        % update stresses
        obj.EnrichGaussDict{1}(ig)=obj.EnrichGaussDict{1}(ig).matsu_enriched(us,ue,ps,pe);
        % calculate the maximum principal stress and its orientation
        % [psmax,gtheta]=obj.EnrichGaussDict{id}.psu;
    end
    for ig=1:length(obj.LineGaussDict{1})
        % update stresses at the discontinuity
        obj.LineGaussDict{1}(ig)=obj.LineGaussDict{1}(ig).matsu_enriched(us,ue,ps,pe);
        % calculate the maximum principal stress and its orientation
        % [psmax,gtheta]=obj.EnrichGaussDict{id}.psu;
    end
% end
% Store the results to corresponding nodes
% May use the nearest three gaussian points instead of the nearest two to
% approximate the stress value 10/06/20. After Khoei 2015 (section 5.6.5.2)
for in=1:obj.NoNodes
    dist=sqrt((xg-obj.X(in)).^2+(yg-obj.Y(in)).^2);
    [~,I]=sort(dist);   % assended order
    stressp=(obj.EnrichGaussDict{1}(I(1)).Stressp+...
        obj.EnrichGaussDict{1}(I(2)).Stressp+...
        obj.EnrichGaussDict{1}(I(3)).Stressp)/3;
    stress=(obj.EnrichGaussDict{1}(I(1)).Stress+...
        obj.EnrichGaussDict{1}(I(2)).Stress+...
        obj.EnrichGaussDict{1}(I(3)).Stress)/3;
    % column stress vector: [rep;sx;sy;sxy;sz] 
    obj.NodDict(in).Stressp=obj.NodDict(in).Stressp+[1;stressp];
    obj.NodDict(in).Stress=obj.NodDict(in).Stress+[1;stress];
end
end

