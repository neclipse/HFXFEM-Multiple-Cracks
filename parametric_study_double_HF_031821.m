% retrieve mesh for all
mesh=postdict(1).EnrichItems{1}.Mygeo.Mesh;
mesh.plotmesh;
hold on;
% specify the ind
ind=56;
crackdict=[];
% retrieve the crackdict
for ienr=1:length(postdict(1).EnrichItems)
    crackgeo=postdict(ind).EnrichItems{ienr}.Mygeo;
    crackdict=[crackdict,crackgeo];
    crackgeo.plotme; % deform, crack, node, phi,ux,uy,scale
end
hold off;
