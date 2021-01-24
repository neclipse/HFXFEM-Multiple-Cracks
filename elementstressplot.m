%% stress, pressure and disp in one encriched elemeent during debugging
allgauss=obj.INTELEM(1,4).EnrichGauss;
% allgauss=obj.EnrichItems(1).INTELEM(1,2).EnrichGauss;
% allgauss=obj.Mytips.NextElem.GaussPntDictM;
% elems=mesh.findelems([6.35,30],elemdict);
% elem=elemdict(elems(1));
% elem.calstress;
% allgauss=elem.GaussPntDictM;
% allgauss=obj.Elemahead.EnrichGauss;
x=[allgauss.X];
y=[allgauss.Y];
p=[allgauss.P];
u=[allgauss.U];
% phi=zeros(size(allgauss));
% for i=1:length(allgauss)
%     phi(i)=allgauss(i).Enf{2}.Phi;
% end
ux=u(1,:);
uy=u(2,:);
spn=[allgauss.Stressp];
figure
plot3(x,y,spn(1,:),'o')
figure
plot3(x,y,spn(2,:),'o')
% figure
% plot3(x,y,phi,'o')
figure
plot3(x,y,p,'*')
figure
plot3(x,y,ux,'s')
figure
plot3(x,y,uy,'s')
%% Common postprocessing
% obj.LeakoffVolume/obj.CrackVolume
% encrack.LeakoffVolume/encrack.CrackVolume
% the second parameter is mesh structure, "3" means trimesh. In fact, I can
% put any number, as long as its size is not 2, which is for quadmesh.
% Step1.Postprocess(6).plotsurf('SPN',3,0,2);
%Step1.Postprocess(end).plotsurf('SPN',3,0,2);
% Step1.Postprocess(end).plotsurf('P',3,0);
% obj.Postprocess(4).plotsurf('P',3,0);
% obj.Postprocess(4).plotsurf('SPN',3,0,2);
% obj.showme(1,'Pfrack', 'second', 'Aperture');
% obj.showme(1,'Pfrack','second','CTraction','component',2);
% obj.showme(1,'Aperture','second','CTraction','component',2);
postprocess.EnrichItems{1}.showme(1,'Aperture','second','CTraction','component',2);
obj.EnrichItems(1).showme(1,'Aperture')
obj.EnrichItems(2).showme(3,'Aperture')

%% Plot crack shape history at give time
Step1.snapshot('crackgeo',47,'timeincs',30,'basecoordinate',30);
% Step1.snapshot('crackgeo',[1,10,20,40]);
timelist=[0,5,10,15,20,25,30]; % time steps to report
timeincs=[Step1.Postprocess.Inc];
incmat=repmat(timeincs,length(timelist),1);
timediff=incmat-transpose(timelist);
[~,inclist]=min(abs(timediff),[],2);
Step1.snapshot('crackgeo',inclist,'timeincs',timelist);

%% element plot for arbitrary element

elems=obj.Preprocess.Mesher.findelems([0.08,30.04],obj.ElemDict);
elem=obj.ElemDict(elems(1));
elem.calstress(false);
allgauss=elem.GaussPntDictM;
% allgauss=obj.Elemahead.GaussPntDictM;
x=[allgauss.X];
y=[allgauss.Y];
p=[allgauss.P];
u=[allgauss.U];
ux=u(1,:);
uy=u(2,:);
% spn=[allgauss.Stressp];
% figure
% plot3(x,y,spn(2,:),'o')
figure
plot3(x,y,p,'*')
% figure
% plot3(x,y,uy,'s')
% figure
% plot3(x,y,ux,'s')
%%
elem=obj.Elemdict(9253);
elem.calstress;
allgauss=elem.GaussPntDictM;
% allgauss=obj.Elemahead.GaussPntDictM;
x=[allgauss.X];
y=[allgauss.Y];
p=[allgauss.P];
uy=[allgauss.Uy];
spn=[allgauss.Stressp];
figure
plot3(x,y,spn(2,:),'o')
figure
plot3(x,y,p,'*')
%% path plot

