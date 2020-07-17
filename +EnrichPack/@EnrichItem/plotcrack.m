function  plotcrack( obj,varargin )
% Author: Chang Huang. 04152019
%PLOTCRACK Summary of this function goes here
%   Detailed explanation goes here
p=inputParser;
defaultdeform=true;
enrichnode=false;
defaultUX=zeros(size(obj.Mesh.VX));
checkUX= @(x) length(x)==length(obj.Mesh.VX) && isnumeric(x);
defaultUY=zeros(size(obj.Mesh.VY));
checkUY= @(x) length(x)==length(obj.Mesh.VY) && isnumeric(x);
defaultscale=[];
checkscale=@(x) isnumeric(x) && x>0;
addOptional(p,'deformflag',defaultdeform);
addOptional(p,'enrichnode',enrichnode);
addOptional(p,'UX',defaultUX,checkUX);
addOptional(p,'UY',defaultUY,checkUY);
addOptional(p,'scale',defaultscale,checkscale);
parse(p,varargin{:});
obj.Mesh.plotmesh(p.Results.deformflag,p.Results.UX,p.Results.UY,p.Results.scale)
hold on
obj.Mygeo.plotme(p.Results.deformflag,0,p.Results.enrichnode,0,p.Results.UX,p.Results.UY,p.Results.scale);
deformedPlus=obj.IntPoints+obj.Uplus*p.Results.scale;
deformedMinus=obj.IntPoints+obj.Uminus*p.Results.scale;
% two lines
% plot(deformedPlus(:,1),deformedPlus(:,2),'w-','LineWidth',1);
% plot(deformedMinus(:,1),deformedMinus(:,2),'w-','LineWidth',1);
% fill plot
crackfill=[deformedPlus;flipud(deformedMinus)];
fill(crackfill(:,1),crackfill(:,2),'w')     % the mesh is gray, the crack space is white
hold off
end

