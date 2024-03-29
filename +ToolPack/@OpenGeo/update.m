function  update( obj )
% Not used in the current build 09/16/20
% Updated calculation of signed distance of nodes within the narrow
% band arond the open geometry
%   1. This method will first use enough discrete points along the crack
%   and the search radius to find all possible interacting elements and
%   nodes
%   2. Call distance2curve function to calculate the shortest distances and
%   the closest points on the curve (x*,y*). 
%   2. Alternatively, I can write the function similar to Pais' code
%   3. Use (x*,y*) to locate the interacted elements using isinside
%   function in FEPack.Elem2d
%   4. Find the nearest two points to (x*,y*) on the curve (x1,y1) and (x2,y2), 
%   sign=(y2-y1)*(x*-x))+(x1-x2)*(y*-y), psi=sign*dist

%  1'. Rediscretize the obj.Segments using proper number of points
    obj.callength;
    obj.calelen;
	maxelen=obj.Maxelength;
	searchr=2.5*maxelen;					% Search radius of possible nodes
	dislen=min(searchr,obj.Length);		% obj.Length must be longer than two element lengths
	np=ceil(obj.Length/dislen)+1;		% number of points to search for the narrow band
    if np>size(obj.Segments,1)
        obj.discretize(np);					% Rediscretize the segments using the desired np
    end
    np=size(obj.Segments,1);
%  1'' Search for the possible  using a loop over all newly generated segments points
	possiblenodes=zeros(100*np,1);
    iloc1=1;
    for ip=1:np
        x=obj.Segments(ip,2);
        y=obj.Segments(ip,3);
        X=obj.Mesh.VX;
        Y=obj.Mesh.VY;
        xdis=X-x;
        ydis=Y-y;
        dis=sqrt(xdis.^2+ydis.^2);
        foundnodes=find(dis<searchr);
        iloc2=iloc1+length(foundnodes)-1;
        possiblenodes(iloc1:iloc2)=foundnodes;
        iloc1=iloc2+1;
    end
    possiblenodes=possiblenodes(1:iloc2);
	nodespool=unique(possiblenodes);	% Eradicate the duplicated nodes in the possiblenodes
%	2. Call distance2curve function to calculate the signed distances and
%   the closest points on the curve (x*,y*)
	plist=obj.Mesh.p(nodespool,:); 		% Retrieve the coordinates of selected nodes
	curvexy=obj.Segments(:,2:3);
	% if obj.Linetype==1
		% interpmethod='Linear';
	% elseif obj.Linetype==2
		% interpmethod='Linear';
	% end
	[xy,distance,~,minind] = distance2curve(curvexy,plist,'linear');
    % xy - the closest point (x*,y*) on the curve to plist
    % distance is the shortest distance from plist to the curve
    % minind - the index of the segment where xy belongs, note minind may
    % be larger than length(obj.Segments) at the last segment, adjustment is
    % needed for the last segment
    L=minind==np;
    if any(L)
        minind(L)=np-1;
    end
    % calculate signed distances, Phi
    x1=obj.Segments(minind,2);
    y1=obj.Segments(minind,3);
    x2=obj.Segments(minind+1,2);
    y2=obj.Segments(minind+1,3);
    xs=xy(:,1);
    ys=xy(:,2);
    x=plist(:,1);
    y=plist(:,2);
    % sign(normal_vec*vec(x*-x))
    signd=sign((y2-y1).*(xs-x)+(x1-x2).*(ys-y));
    Phipool=distance.*signd;                                                % signed distance function of all nodes in the pool
    L2=Phipool==0;
    if any(L2)
       Phipool(L2)=1e-12;                                                   % Kick nodes off the crack line 
    end
%   3. Use (x*,y*) to locate the interacted elements using isinside
    %   function in FEPack.Elem2d
    bodyelems=obj.findelems(xy);                                            
%   4. Find the RTipelements and Bodyelements
    rtipelems=obj.findelems(obj.Tips(obj.Rtips,:));
    obj.Rtipelements=rtipelems;
    % find the Bodyelements by erradicating the tip elements
    obj.Bodyelements=setdiff(bodyelems,rtipelems);
%   5. Assign Phi to directly interacted nodes
    numnode=8*length(bodyelems);
    nodes=zeros(numnode,1);                                                 % preallocate the array for the directly interacted nodes
    iloc1=1;
    for ielem=1:length(bodyelems)
        nodlist=obj.Elemdict(bodyelems(ielem)).NodList;
        iloc2=iloc1+length(nodlist)-1;
        nodes(iloc1:iloc2)=nodlist;
        iloc1=iloc2+1;
    end
    nodes=nodes(1:iloc2);
    Nodes=unique(nodes);                                                    % indices of the directly interacted nodes
    Phi=zeros(length(Nodes),2);                                             % first col has the node ind, and second col has the Phi
    Phi(:,1)=Nodes;
    for i=1:length(Nodes)
        Phi(i,2)=Phipool(nodespool==Nodes(i));
    end
    obj.Nodes=Nodes;
    obj.Phi=Phi;                                                            % The signed distance function, ie., level set function phi
%   5' Find tip nodes    
    numnode=8*length(rtipelems);
    nodes=zeros(numnode,1);
    iloc1=1;
    for ielem=1:length(rtipelems)
        nodlist=obj.Elemdict(rtipelems(ielem)).NodList;
        iloc2=iloc1+length(nodlist)-1;
        nodes(iloc1:iloc2)=nodlist;
        iloc1=iloc2+1;
    end
    nodes=nodes(1:iloc2);
    obj.Rtipnodes=unique(nodes);
    obj.Bodynodes=setdiff(obj.Nodes,obj.Rtipnodes);
%     obj.findblending(2,1);
end

