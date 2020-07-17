function h=plotmesh_old( obj,varargin )
%PLOTMESH, METHOD OF QUADMESHER
%reshape the vectors to be compatible with the requirement of the built-in 'surf'
% if ~isempty(varargin)
%     deformflag=varargin{1};
%     UX=varargin{2};
%     UY=varargin{3};
% else
%     deformflag=false;
% end
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

%%
meshstructure=[obj.tannodes,obj.radnodes];
if ~p.Results.deformflag
    X=transpose(reshape(obj.VX,meshstructure));
    Y=transpose(reshape(obj.VY,meshstructure));
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
    Xd=obj.VX+p.Results.UX*scale;
    Yd=obj.VY+p.Results.UY*scale;
    X=transpose(reshape(Xd,meshstructure));
    Y=transpose(reshape(Yd,meshstructure));
    titlename='Mesh plot of the deformed domain';
end
    var=transpose(reshape(obj.Nanfield,meshstructure));
    figure('Name',titlename,'NumberTitle','off');
    h=surf(X,Y,var);
    axis equal % equal unit length along the x- and y- directions.
    % colormap (jet)
    % colorbar
    colormap gray
    title(titlename)
    xlabel('x')
    ylabel('y')
    view(0,90)
end


