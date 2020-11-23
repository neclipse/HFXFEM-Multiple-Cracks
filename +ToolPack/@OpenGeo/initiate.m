function  initiate( obj,varargin )
% Initial calculation of signed distance of nodes within the narrow
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

    searchlevel=3;
% if nargin<3
%     option=true;        % by default, apply Toedge
% else
%     option=varargin{2};
% end

%   1. use obj.Toedge_1 to extend the tips to element edge if needed
    rtipelems1=obj.Toedge1;
    %  1'. Rediscretize the obj.Segments using proper number of points
    obj.callength;
    obj.calelen;
	maxelen=obj.Maxelength;
	searchr=searchlevel*maxelen;		% Search radius of possible nodes
	dislen=min(searchr,obj.Length);		% obj.Length must be longer than two element lengths
	np=ceil(obj.Length/dislen)+1;		% number of points to search for the narrow band
    if np>size(obj.Segments,1)
        obj.discretize(np);					% Rediscretize the segments using the desired np
    end
    np=size(obj.Segments,1);
%  1'' Search for the possible nodes using a loop over all newly generated segments points
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
	[xy,~,~,~] = distance2curve(curvexy,plist,'linear');
    % (distance) is the shortest distance from plist to the curve (excluding
    % the extension of the curve, consequently many xy will be the endpoint 
    % if plist cannot find its normal projection within the
    % curve).obsolete.
    % see the following change about extendedcurve.
    
    % It is later on deemed important to find normal projection and real
    % closest distance. So we extend the curve on 06292019, the change
    % extends the definition of y2.
    r=2*obj.Maxelength;
    extendedtip1x=obj.Tips(1,1)+r*cos(obj.Omegas(1));
    extendedtip1y=obj.Tips(1,2)+r*sin(obj.Omegas(1));
    extendedtip2x=obj.Tips(2,1)+r*cos(obj.Omegas(2));
    extendedtip2y=obj.Tips(2,2)+r*sin(obj.Omegas(2));
    extendedcurvexy=[extendedtip2x,extendedtip2y;curvexy;extendedtip1x,extendedtip1y];
    [xytemp,distance,~,minind] = distance2curve(extendedcurvexy,plist,'linear');
    % xy - the closest point (x*,y*) on the curve to plist
    % distance is the shortest distance from plist to the curve (including 
    % the extension)
    % minind - the index of the segment where xy belongs, note minind 
    % is length(obj.Segments)+1 at the last segment, adjustment is
    % needed for the last segment
    if length(distance(distance<1e-2*obj.Minelength))>2
        warning('mesh coincides with No. %d crack, please revise mesh',obj.Id);
    end
    % calculate signed distances, Phi
%     % Special care to the last point of the (extended)crack 06292019 change
%     % np to np+2 because extended curve has two more points.
    minind(minind==np+2)=np+1;
    x1=extendedcurvexy(minind,1);
    y1=extendedcurvexy(minind,2);
    x2=extendedcurvexy(minind+1,1);
    y2=extendedcurvexy(minind+1,2);
    xs=xytemp(:,1);
    ys=xytemp(:,2);
    x=plist(:,1);
    y=plist(:,2);
    % sign(normal_vec*vec(x-x*))
    signd=sign((y1-y2).*(x-xs)+(x2-x1).*(y-ys));
    Phipool=distance.*signd;                                                % signed distance function of all nodes in the pool
    L2=(Phipool==0);
    if any(L2)
       Phipool(L2)=1e-8;                                                    % Kick nodes off the crack line 
    end
    obj.Nodespool=nodespool;
    obj.Phipool=Phipool;
%   3. Use (x*,y*) to locate the interacted elements using isinside
    %   function in FEPack.Elem2d
    xy=uniquetol(xy,'ByRows',true);                                                   % added on 1/31 to avoid extra work on the same point
    interactedelems=obj.findelems(xy,'initial');                            % changed from 'inside' to 'initial' to have count more elements
    obj.Intelements=interactedelems;

%   4. Delete bodyelements with tiny impact area, then goodelements left
    % Find the intersections
    % MISTAKE IN TESTING IF A INTERACTED ELEMENT IS GOOD: CANNOT TELL WHEN THE CRACK IS
    % JUST TOUCHING AN ELEMENT (01/31/2019), FIXED on 02/01/2019 by adding
    % constraints in findelemes
    goodlist=true(size(interactedelems));
    for ielem=1:length(interactedelems)
        elem=obj.Elemdict(interactedelems(ielem));
        goodlist(ielem)= intersection(obj, elem);
%   The following operations are correct but should not be function of
%   mygeo. Therefore, it is preferred to be done in
%   EnrCrackBody.initial_enrich although it will require a second call of
%   mygeo.intersection method. 11/20/2020
%         if goodlist(ielem)
%             elem.setenrich(obj.Id); 
%             ind=elem.Enrich==obj.Id;
%             elem.LocalInt{ind}=localpnts;
%             elem.GlobalInt{ind}=pnts;
%         end
    end
    % BUG: the following line will add back the elements filtered out by
    % intersection method. Need to be commented. 07/07/2020.
    %goodlist=true(size(interactedelems));
    goodelements=interactedelems(goodlist);
    obj.Intelements=goodelements;
%   5. Assign Phi to directly interacted nodes
    numnode=8*length(goodelements);
    nodes=zeros(numnode,1);                                                 % preallocate the array for the directly interacted nodes
    iloc1=1;
    for ielem=1:length(goodelements)
        nodlist=obj.Elemdict(goodelements(ielem)).NodList;
        iloc2=iloc1+length(nodlist)-1;
        nodes(iloc1:iloc2)=nodlist;
        iloc1=iloc2+1;
    end
    nodes=nodes(1:iloc2);
    Nodes=unique(nodes);                                                    % indices of the directly interacted nodes
    Phi=zeros(length(Nodes),2);                                             % first col has the node ind, and second col has the Phi
    Phi(:,1)=Nodes;
    if ~isempty(setdiff(Nodes,nodespool))
       error('The specified searchlevel may be too small, please increase it to 3.5 or larger') 
    end
    for i=1:length(Nodes)
        % nodespool and Phipool are of the same order
        Phi(i,2)=Phipool(nodespool==Nodes(i));
    end
    obj.Nodes=Nodes;
    obj.Phi=Phi;                                                            % The signed distance function, ie., level set function phi  
%   6. use obj.Toedge_2 to find the nodes on the tipped element edges
    rtipelems2=obj.Toedge2;
    rtipelems=[rtipelems1;rtipelems2];
    rtipelems(rtipelems==0)=[];
    rtipelems=unique(rtipelems);
%     obj.Phipool=[nodespool,Phipool];
%   7. Find tip nodes (all nodes in the tip elements)
    numnode=8*length(rtipelems);
    if numnode>0
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
    end
    obj.Bodynodes=setdiff(obj.Nodes,obj.Rtipnodes);
    obj.Rtipelements=rtipelems;
    obj.Bodyelements=setdiff(goodelements,rtipelems);
    %   8. Find the injection point surrounding elements (blending elements)
    if ~isempty(varargin)
        injectionpoint=varargin{1};
        obj.findblending(1,injectionpoint,5);
    end
end




































