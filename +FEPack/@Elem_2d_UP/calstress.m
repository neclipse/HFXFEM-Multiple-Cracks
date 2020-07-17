function calstress( obj)
%Calstress: function of elem2d, will be called in domain.storage
%Extrapolate stress from conventional gaussian points to the nodes
% Also obtain the stress at the elemental centroid
xi=[-1;1;1;-1;0];                   % Four nodes and the centroid in a column vector
eta=[-1;-1;1;1;0];
[stressp,stress]=obj.extraplstress(xi,eta); 
% stress at the element centroid
obj.Stress=stress(end,:)';
obj.Stressp=stressp(end,:)';
% Store the results to corresponding nodes 
for in=1:obj.NoNodes
    % column stress vector: [rep,sx,sy,sxy,sz] 
    obj.NodDict(in).Stressp=obj.NodDict(in).Stressp+[1;stressp(in,:)'];
    obj.NodDict(in).Stress=obj.NodDict(in).Stress+[1;stress(in,:)'];
end
end

