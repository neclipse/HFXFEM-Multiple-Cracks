function [fe,locarray]=cal_qextenr(obj,q,varargin)
%CAL_QEXTENR to incorporate crack mouth inflow. 03172019.
% Since 04292019, this function actually calculate both the standard and
% the enriched parts of qext, so that one do not need specify Qtable_node
% for preprocessor.

% 12/01/19Note the element edge is treated as an external edge, we do not
% need to do line integral again over the neighboring element. This is
% equivalent to divide the q by half and integrate twice along the same
% edge.

% use the input flow at the crack mouth and local nodes index to calculate
% the injection volume to the enriched dof
%   1. locate crack mouth using the obj.Mygeo (the tip that is not rtip)
%   or locate the element using the varargin{1}, which is a coordinate of
%   the injection point.
%   2. Go to the enriched element and the specified edge, and then
%   specify the gauss nodes along the edge
%   3. Use the obj. Myenf to calculate the ridge function value of the
%   gauss point
%   4. use the line integral equation to calculate the fe
%   5. Also return the locarray_enr of the specified nodes
feenr=zeros(2,1);
festd=zeros(2,1);
% agl=assembleglobalinputs();
% skin=agl.skin;
% 1. Find mouthelemid and the two closest nodes from the injection point
if isempty(varargin)
    tipid=setdiff([1,2],obj.Mygeo.Rtips);
    point=obj.Mygeo.Tips(tipid,:); %[x,y]
else
    point=varargin{1}; %[x,y]
end
interactedelems=obj.Mygeo.findelems(point,'in_edge');
mouthelemid=intersect(interactedelems,obj.Mygeo.Intelements);
mouthelem=obj.Elemdict(mouthelemid);
dist=sqrt((mouthelem.X-point(1)).^2+(mouthelem.Y-point(2)).^2);
[~,I]=sort(dist,'ascend');
nodes=I(1:2);
% 2. Find EnFRidge and get the nodes_phi
for ienf=1:length(obj.Myenfs)
    if isa(obj.Myenfs{ienf},'EnrichPack.EnFRidge')
        enfridge=obj.Myenfs{ienf};
        continue;
    end
end
allnodes=mouthelem.NodList;
nodespool=enfridge.Lsv(:,1);                  % all nodes in the interacted elems
phipool=enfridge.Lsv(:,2);
[Lia,Lcob]=ismember(allnodes,nodespool);    % find the relative index of nodes in the nodespool
if ~all(Lia)
    error('The requested nodes are not stored in the nodespool');
end
nodes_phi=phipool(Lcob);
% 2. To the edge gauss points
if mouthelem.NoNodes==4
parentlocals=[-1,-1;1,-1;1,1;-1,1];   % local coordinates [ksi,eta] of four nodes
end
% X=mouthelem.X(nodes);
% Y=mouthelem.Y(nodes);
locals=parentlocals(nodes,:);         % local coordinates of the given two nodes
p=1;
[xx,ww]=GaussQuad(p); 
Nlinematrix=[1-xx,1+xx]/2;
gp=Nlinematrix*locals;                     % two dimensional gp
gw=ww;
gaussdictm(p)=FEPack.GaussPnt_LE_UP();
% 3. Calculate the ridge function on the selected 2D gp
% 4. Do the line integral using 1D [xx,ww], which is essentially the same
% as 2D gp but no need to drop unnecessary zeros.
% if mode==1
%     % 5. Find the crack opening for injection point, 0429 thinking how to 
%     % apply concentrated flow to the crack.
%     ws=obj.Aperture;
%     intpoints=obj.IntPoints;
%     tip=intpoints(1,:);
%     temp=intpoints-repmat(tip,size(intpoints,1),1);
%     dist=sqrt(temp(:,1).^2+temp(:,2).^2);
%     distq=sqrt((point(1)-tip(1))^2+(point(2)-tip(2))^2);
%     wq=interp1(dist,ws,distq);
%     if wq==0
%        wq=1e-6; % to initiate flow through the crack 
%     end
%     q=q*wq;
% end
gbinp=obj.Elemdict(1).GaussPntDictM(1).GBINP;
for i=1:p
gaussdictm(i)=FEPack.GaussPnt_LE_UP(gp(i,1),gp(i,2),gw(i),mouthelem.NoNodes,gbinp);
gaussdictm(i)=gaussdictm(i).preparing(mouthelem.X,mouthelem.Y);
Np=gaussdictm(i).Np;
ridge=enfridge.calculate(nodes_phi,Np);
xi=xx(i);                            % one dimensional gp
[Nline,~]=lineshape(mouthelem.NoNodes,xi);
% [Nline,Nline_xi]=lineshape(mouthelem.NoNodes,xi);
% x_xi=Nline_xi*X;
% y_xi=Nline_xi*Y;
% From 09152019, the injection rate becomes injection volume, so ds is not
% needed. Otherwise, the actual injeciton volume will be dependent on the 
% element size.
% ds=sqrt(x_xi^2+y_xi^2);
ds=1/2;
Nlinep=[Nline(1,1),Nline(1,3)];
feenr=feenr+ww(i)*Nlinep'*ridge*q*ds;
% For later on,  qinjection is calculated in this single function, 0429
festd=festd+ww(i)*Nlinep'*q*ds;
fetot=festd+feenr;
fe=[festd;feenr];
% fe=fetot;
end
obj.InjectionFlux=-fetot;        % positive is injection.
for in=1:length(nodes)
    mouthelem.NodDict(nodes(in)).CInjection=-fetot(in);
end
% 5. Find the locarray_enr and locarray_std
% ERROR found on 0423  FIXED on 0424, do not use LocarraryPEnr
locarray_enr=mouthelem.JacobianMatDict(obj.Id).LocarrayEnr(3*nodes);
locarray_std=mouthelem.Locarray(3*nodes);
locarray=[locarray_std,locarray_enr];
% locarray=locarray_enr;
end

