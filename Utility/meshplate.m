function [ node,element,wnode,hnode ] = meshplate( lw,lh,lc,xstep,ystep,varargin )
%MESHPLATE mesh preparation for nodes info and element info,
% by author: Chang Huang at LSU (huangchang73@gmail.com)
% Quadrilateral plate with a center crack
% only needs to generate xv and yv
% refer to meshgrid, and connectivity_nonuniform
% the connection can be defined easily using the global index of nodes

% inputs:
% lw: width
% lh: height
% lc: crack length
% xstep: total element number in one row along x-direction
% ystep: total element number in one column along y-direction
% bias: The bias ratio from reference element length to largest element
% length, which equals to 1, if not stated.

%% generate the structured coordinates array, xgv and ygv
wnode=xstep+1;                      % xstep is the section number
hnode=ystep+1;
if ~isempty(varargin)                   % Nonuniform
    bias=varargin{1};                  % Lxelem/lxelem;
    xcrackstep=varargin{2};              % element number in the dense zone around the crack
    ycrackstep=varargin{3};
    xratio=varargin{4};                  % the ratio of the length of the refined mesh zone to the crack length
    yratio=varargin{5};
    % the center densed zone
    xgv1=linspace(-xratio*lc/2,xratio*lc/2,xcrackstep+1);
    ygv1=linspace(-yratio*lc/2,yratio*lc/2,ycrackstep+1);
    % calculate the loose zone
    swidth=(lw-xratio*lc)/2;                 % length of the loose part
    sheight=(lh-yratio*lc)/2;
    sxstep=(xstep-xcrackstep)/2;         % section number in the half loose part
    systep=(ystep-ycrackstep)/2;
    % scale factor
    qx=bias^(1/(sxstep-1));
    qy=bias^(1/(systep-1));
    % calculate the first element length
    xl1=swidth*(1-qx)/(1-qx^sxstep);    % from the sum of geometric array
    yl1=sheight*(1-qy)/(1-qy^systep);
    % calculate the coordinates array 
    % right part
    xgv2=zeros(1,sxstep);
    ygv2=zeros(1,systep);
    xcoord=xratio*lc/2;
    ycoord=yratio*lc/2;
    for ind=1:sxstep
        seg=xl1*qx^(ind-1);
        xcoord=xcoord+seg;
        xgv2(ind)=xcoord;
    end
    for ind=1:systep
        seg=yl1*qy^(ind-1);
        ycoord=ycoord+seg;
        ygv2(ind)=ycoord;
    end
    % left part and right part is symmetric about the center axis 
    xgv3=fliplr(-xgv2);
    ygv3=fliplr(-ygv2);
    % assemble the coordinates array
    xgv=[xgv3,xgv1,xgv2];
    ygv=[ygv3,ygv1,ygv2];
else                                    % Uniform
    xgv=linspace(-lw/2,lw/2,wnode);
    ygv=linspace(-lh/2,lh/2,hnode);    
end
%% generate the 2d mesh points, p or node
[X,Y]=meshgrid(xgv,ygv);
NNode  = size(X,1)*size(X,2);
nXElem = length(xgv)-1;
nYElem = length(ygv)-1;
XYZ   = zeros(NNode,2);
NN = zeros(wnode,hnode);
nNode=1;
for iXNode = 1:(nXElem+1)
    for iYNode = 1:(nYElem+1)
        NN(iXNode,iYNode) = nNode;
        XYZ(nNode,:) = [ X(iYNode,iXNode), Y(iYNode,iXNode)];
        nNode = nNode+1;
    end
end    
%% generate the global connection table, t or EToV or element
nElem = 1;
CONNEC = zeros(nXElem*nYElem,4);
for iXElem = 1:nXElem
    for iYElem = 1:nYElem
        N1 = NN(iXElem,iYElem);
        N2 = NN(iXElem+1,iYElem);
        N3 = NN(iXElem+1,iYElem+1);
        N4 = NN(iXElem,iYElem+1);
        CONNEC(nElem,1:4) = [N1 N2 N3 N4];
        nElem = nElem+1;
    end
end
node=XYZ;
element=CONNEC;

