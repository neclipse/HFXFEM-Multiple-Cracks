function initiate(obj,varargin)
% global threshold;
% agl=assembleglobalinputs();
agl=obj.Elemdict(1).GaussPntDictM(1).GBINP;
threshold=agl.threshold_formaxps;
obj.Mesh=obj.Mygeo.Mesh;
obj.Xct=obj.Mygeo.Tips(obj.Itip,1);
obj.Yct=obj.Mygeo.Tips(obj.Itip,2);
obj.Omega=obj.Mygeo.Omegas(obj.Itip);
obj.setinteractedelem(obj.Itip);
obj.setenrichednode;
obj.setradius(2);   % the ratio of search radius over element length
% by default the check and direction are based on maxps
obj.Growcheck=ToolPack.Maxpscheck(threshold,'tolerance',0.05,'tolerance2',0.15);
% obj.Growdirection=ToolPack.Maxpsdirection(obj.Omega);
end