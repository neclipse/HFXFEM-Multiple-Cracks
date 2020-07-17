function  plotcrack( obj,mesh,varargin )
%PLOTCRACK is a function of postprocessor, similar to enrichitem.plotcrack
%   Detailed explanation goes here

% Author: Chang Huang. 07152019

p=inputParser;
defaultdeform=true;
enrichnode=false;
defmeshflag=false;
defaultscale=1;
checkscale=@(x) isnumeric(x) && x>0;
addOptional(p,'meshflag',defmeshflag);
addOptional(p,'deformflag',defaultdeform);
addOptional(p,'scale',defaultscale,checkscale);
addOptional(p,'enrichnode',enrichnode);


parse(p,varargin{:});
if p.Results.meshflag
    mesh.plotmesh(p.Results.deformflag,obj.UX,obj.UY,p.Results.scale)
    hold on
end
if p.Results.deformflag
    scale=p.Results.scale;
else
    scale=defaultscale;
end
% obj.Mygeo.plotme(p.Results.deformflag,0,p.Results.enrichnode,0,obj.UX,obj.UY,p.Results.scale);
for i=1:length(obj.EnrichItems)
deformedPlus=obj.EnrichItems{i}.IntPoints+obj.EnrichItems{i}.Uplus*scale;
deformedMinus=obj.EnrichItems{i}.IntPoints+obj.EnrichItems{i}.Uminus*scale;
% two lines
% plot(deformedPlus(:,1),deformedPlus(:,2),'w-','LineWidth',1);
% plot(deformedMinus(:,1),deformedMinus(:,2),'w-','LineWidth',1);
% fill plot
crackfill=[deformedPlus;flipud(deformedMinus)];
fill(crackfill(:,1),crackfill(:,2),'w','LineWidth',1)     % the mesh is gray, the crack space is white
hold on
end
hold off
end





