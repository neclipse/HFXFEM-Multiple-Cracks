function calstress_enriched(obj,storage)
%Calstress_enriched: function of elem2d, will be called in domain.storage
%Extrapolate/Average stress from subtriangles gaussian points to the nodes
%1. Loop over each node to find the closest three gaussian points
%2. Obtain the average stresses as the contribution of stresses from 
% this element to the node, only when storage==true
% obtain the coordinates of the enriched gaussian points
xg=[obj.EnrichGauss.X];
yg=[obj.EnrichGauss.Y];
%% Calculate the stress at the gaussian points, they are stored for this
% increment. That's why Maxpscheck do not need to call this method again.
us=obj.Un2i;      % elemental iterative total displacement vector, u_sup(n+1)_sub(i+1)
ps=obj.Pn2i;      % elemental iterative total pore pressure vector, p_sup(n+1)_sub(i+1)
Jacob=obj.JacobianMat;
ue=Jacob.Un2iEnr;
pe=Jacob.Pn2iEnr;
for ig=1:length(obj.EnrichGauss)
    % update stresses
    obj.EnrichGauss(ig)=obj.EnrichGauss(ig).matsu_enriched(us,ue,ps,pe);
    % calculate the maximum principal stress and its orientation
    % [psmax,gtheta]=obj.EnrichGaussDict{id}.psu;
end
% The stress along the linegauss may not be accurate as this is exactly on
% the discontinuity, the stress changes drastically. 
%for ienr=1:obj.EnrichNum
%     for ig=1:length(obj.LineGaussDict{ienr})
%         % update stresses at the discontinuity
%         obj.LineGaussDict{ienr}(ig)=obj.LineGaussDict{ienr}(ig).matsu_enriched(us,ue,ps,pe);
%         % calculate the maximum principal stress and its orientation
%         % [psmax,gtheta]=obj.EnrichGaussDict{id}.psu;
%     end
% end
%% Approximate the stress at the element centroid or potential intersection
% Skipped for now. May use nonlocal tip stress for smoothier prediction.
% 12/22/2020

%% Store the results to corresponding nodes
% May use the nearest three gaussian points instead of the nearest two to
% approximate the stress value 10/06/20. After Khoei 2015 (section 5.6.5.2)
% Modified on 12/22/2020
if storage 
    for in=1:obj.NoNodes
        dist=sqrt((xg-obj.X(in)).^2+(yg-obj.Y(in)).^2);
        [~,I]=sort(dist);   % assended order
        stressp=(obj.EnrichGauss(I(1)).Stressp+...
            obj.EnrichGauss(I(2)).Stressp+...
            obj.EnrichGauss(I(3)).Stressp)/3;
        stress=(obj.EnrichGauss(I(1)).Stress+...
            obj.EnrichGauss(I(2)).Stress+...
            obj.EnrichGauss(I(3)).Stress)/3;
        % column stress vector: [rep;sx;sy;sxy;sz]
        obj.NodDict(in).Stressp=obj.NodDict(in).Stressp+[1;stressp];
        obj.NodDict(in).Stress=obj.NodDict(in).Stress+[1;stress];
    end
end
end

