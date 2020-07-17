function update(obj)
obj.Xct=obj.Mygeo.Tips(obj.Itip,1);
obj.Yct=obj.Mygeo.Tips(obj.Itip,2);
obj.Omega=obj.Mygeo.Omegas(obj.Itip);
obj.setinteractedelem(obj.Itip);
obj.setenrichednode;
obj.setradius(2);
obj.checkactive; % not finished 06212019
end