function  gaussdictm=crtgpd( obj, varargin )
%CRTGPD  create the gaussian point dictionary, called by the constructor of
%Elem_2d_EP
% Inputs explanation
% obj: the domain object
% p:   the order of the Gaussian quadrature (p<=12)
p=2;
%% Specify the inputs
if nargin==2
    p=varargin{1}; % the third input is the order of the integreted polynomial
end
%% Specify the required number of gaussian points
nnodes=4;
[Q,W]=gausstable(p,'QUAD');
map=[1,3,4,2];
gauss=[Q,W];
gauss=gauss(map,:);
numgauss=size(gauss,1); % number of gauss points
%% Create proper gaussian point dictionary corresponding to material type
if obj.MatType==1
    gaussdictm(1,numgauss)=FEPack.GaussPnt_LE(); %generate an void object array
    for igauss=1:numgauss
        gaussdictm(igauss)=FEPack.GaussPnt_LE(gauss(igauss,1),gauss(igauss,2),gauss(igauss,3),nnodes,obj.GBINP);
    end
elseif obj.MatType==2
    gaussdictm(1,numgauss)=FEPack.GaussPnt_DP(); %generate an void object array
    for igauss=1:numgauss
        gaussdictm(igauss)=FEPack.GaussPnt_DP(gauss(igauss,1),gauss(igauss,2),gauss(igauss,3),nnodes,obj.GBINP);
    end
elseif obj.MatType==4
    gaussdictm(1,numgauss)=FEPack.GaussPnt_LE_UP(); %generate an void object array
    for igauss=1:numgauss
        gaussdictm(igauss)=FEPack.GaussPnt_LE_UP(gauss(igauss,1),gauss(igauss,2),gauss(igauss,3),nnodes,obj.GBINP);
    end
end
end

