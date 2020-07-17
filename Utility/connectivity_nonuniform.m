% Written By: Matthew Jon Pais, University of Florida (2010)
% Website: http://sites.google.com/site/matthewjpais/Home
% Email: mpais@ufl.edu, matthewjpais@gmail.com

function connectivity_nonuniform
% This function defines connectivity information about the problem.  This
% includes defining node numbers for the problem, defining the XYZ 
% coordinates of the node numbers, creating a connectivity matrix, and
% creating matricies storing the X and Y coordinates for each element.

global CONNEC DOMAIN NODES XYZ

Width   = DOMAIN(1);                                                         % Number of elements in the x-direction
Height  = DOMAIN(2);                                                         % Number of elements in the y-direction
lXElem1 = DOMAIN(3);                                                         % Elemental length in the x-direction
lYElem1 = DOMAIN(4);                                                         % Elemental length in the y-direction
if length(DOMAIN)==8
    cwidth  = DOMAIN(5);                                                     % Length of the center dense zone
    cheight = DOMAIN(6);
    lXElem2 = DOMAIN(7);                                                     % Elemental length in the x-direction in the dense zone
    lYElem2 = DOMAIN(8);                                                     % Elemental length in the y-direction in the dense zone
    swidth  =(Width-cwidth)/2;
    sheight =(Height-cheight)/2; 
    xgv1=linspace(0,swidth,swidth/lXElem1);
    xgv2=linspace(swidth,(swidth+cwidth),cwidth/lXElem2);
    xgv3=linspace((swidth+cwidth),Width,swidth/lXElem1);
    xgv =[xgv1, xgv2(2:end), xgv3(2:end)];
    ygv1=linspace(0,sheight,sheight/lYElem1);
    ygv2=linspace(sheight,(sheight+cheight),cheight/lYElem2);
    ygv3=linspace((sheight+cheight),Height,sheight/lYElem1);
    ygv =[ygv1, ygv2(2:end), ygv3(2:end)];
else
    xgv=0:lXElem1:Width;
    ygv=0:lYElem1:Height;
end
[X,Y]  = meshgrid(xgv,ygv);
NNode  = size(X,1)*size(X,2);
nXElem = length(xgv)-1;
nYElem = length(ygv)-1;
NODES  = zeros(NNode,31);
XYZ    = zeros(NNode,3); 
NN     = zeros(nXElem+1,nYElem+1);
nNode = 1;
for iYNode = 1:(nYElem+1)
    for iXNode = 1:(nXElem+1)
        NN(iXNode,iYNode) = nNode;
        XYZ(nNode,:) = [nNode X(iYNode,iXNode) Y(iYNode,iXNode)];
        NODES(nNode,1) = nNode;
        nNode = nNode+1;
    end
end
% Create a global connectivity matrix (CONNEC) and elemental coordinate matricies (Ex, Ey)
% CONNEC  = [ElementNumber,LocalNode1,LocalNode2,LocalNode3,LocalNode4,InclusionElementNum]
nElem = 1;
CONNEC = zeros(nXElem*nYElem,5);
for iYElem = 1:nYElem
    for iXElem = 1:nXElem
        N1 = NN(iXElem,iYElem);
        N2 = NN(iXElem+1,iYElem);
        N3 = NN(iXElem+1,iYElem+1);
        N4 = NN(iXElem,iYElem+1);
        CONNEC(nElem,1:5) = [nElem N1 N2 N3 N4];
        nElem = nElem+1;
    end
end
DOMAIN(1) = nXElem;                                                         % Number of elements in the x-direction
DOMAIN(2) = nYElem;
end