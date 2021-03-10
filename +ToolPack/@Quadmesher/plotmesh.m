function  plotmesh(obj,varargin)
%PLOTMESH (updated) Using patch function to plot unstructured mesh

%% Inputs Parser and validation added on 04112019.
p=inputParser;
defaultdeform=false;
defaultUX=zeros(size(obj.VX));
checkUX= @(x) length(x)==length(obj.VX) && isnumeric(x);
defaultUY=zeros(size(obj.VY));
checkUY= @(x) length(x)==length(obj.VY) && isnumeric(x);
defaultscale=[];
% checkscale=@(x) isnumeric(x) && x>0;
addOptional(p,'deformflag',defaultdeform);
addOptional(p,'UX',defaultUX,checkUX);
addOptional(p,'UY',defaultUY,checkUY);
addOptional(p,'scale',defaultscale);
addOptional(p,'nodeflag',false);
addOptional(p,'elemflag',false);
parse(p,varargin{:});

%% Plot the mesh using patch function 
% patch will use X, Y to create polygons, each column of X and Y represents
% the coordinates of the vertices of one polygon
X=zeros(obj.Nface,obj.Totelem);
Y=X;
if ~p.Results.deformflag
    Xloc=obj.VX;
    Yloc=obj.VY;
    titlename='Mesh plot of the undeformed domain';
else
    if isempty(p.Results.scale)
        dispmin=min(min([p.Results.UX,p.Results.UY]));
        dispmax=max(max([p.Results.UX,p.Results.UY]));
        maxabsdisp=max(abs(dispmin),abs(dispmax));
        scale=10^(abs(fix(log10(maxabsdisp)))+0.5);      % scaling factor
    else
        scale=p.Results.scale;
    end
    Xloc=obj.VX+p.Results.UX*scale;
    Yloc=obj.VY+p.Results.UY*scale;
    titlename='Mesh plot of the deformed domain';
end
% rearrange into the required structure by patch function
for ielem = 1:obj.Totelem
    nind = obj.EToV(ielem,:);                                                 % Nodes for current element
    X(:,ielem) = Xloc(nind);                                               % Initial x-coordinates of nodes
    Y(:,ielem) = Yloc(nind);                                               % Initial y-coordinates of nodes
end

figure; hold on;
patch(X,Y,[0.8 0.8 0.8])
axis equal
% title(titlename)
xlabel('X'); ylabel('Y'); axis off; axis equal; 

%% Plot the node numbers
if p.Results.nodeflag
    for i = 1:obj.Totnodes
        text(X(i),Y(i),num2str(i));
    end
end

%% Plot the element numbers
if p.Results.elemflag
    for i = 1:obj.Totelem
        nind = obj.EToV(i,:);
        XN_diag = X(nind);
        YN_diag = Y(nind);
        text(mean(XN_diag),mean(YN_diag),num2str(i));
    end
end
end

