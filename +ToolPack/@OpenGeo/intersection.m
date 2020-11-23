function [good,pnts,localpnts] = intersection( obj,elem,varargin )
% intersection method of OpenGeo
% find the intersections of elem
% of the OpenGeo and the element edges
%   Take advantage of the signed distance function
% Outputs:
% 1. good: if the crack has a good impact on the interacted element
% 2. pnts: the global coordinates of the two intersections
% 3. seeds: vertices for delaunay to create triangular subdomains
if ~isempty(varargin)
    nodespool=varargin{1};
    Phipool=varargin{2};
else
    nodespool=obj.Nodespool;
    Phipool=obj.Phipool;
end
corner = [1,2,3,4,1];
node=[-1,-1;1,-1;1,1;-1,1];
good=true;
% seeds=[];                   
pnts=zeros(2,2);            % global coordinates of the intersections
localpnts=zeros(2,2);       % local coordinates of the intersections
% extrapnts=zeros(2,2);       % local coordinates of the two extra points for the seeds
ipnt=1;
pairs=zeros(1,4);           % The nodes on the edges which are intersected by the crack
R=zeros(1,2);
% Find the two intersections and on what edge
for i =1:4
    % global index of node i and i+1
    n1=elem.NodList(corner(i));
    n2=elem.NodList(corner(i+1));
    % global coordinates of node i and node i+1
    xy1=[elem.X(corner(i)),elem.Y(corner(i))];
    xy2=[elem.X(corner(i+1)),elem.Y(corner(i+1))];
    % local coordinates of node i and i+1
    local1=node(corner(i),:);
    local2=node(corner(i+1),:);
    % signed distance of node i, i+1
    phi1=Phipool(nodespool==n1);
    phi2=Phipool(nodespool==n2);
    % find the intersection and its global coordinates
    % May need to be corrected as there may be kinks at the element edges
    if phi1*phi2<0
        pairs(2*ipnt-1:2*ipnt)=[corner(i),corner(i+1)];
        r = phi1/(phi1-phi2);
        R(ipnt)=r;
        % global coordinate of intersections
        pnts(ipnt,:)=(1-r)*xy1+r*xy2;
        % local coordinate of intersections
        localpnts(ipnt,:)=(1-r)*local1+r*local2;
        ipnt=ipnt+1;
    end
end
% REORDER THE INTERSECTIONS TO ALIGN WITH THE DIRECTION OF CRACK 11/05/2018
% The reorderring is to let the crtip.realgrow know which point is the new
% point. 10/30/20. 
% This is also to make sure the elem.linegauss has the right (consistent) crack
% direction 11/04/20.
% Find the nearest segments
 p=xy1;
 curvexy=obj.Segments(:,2:3);
 [~,~,~,minind] = distance2curve(curvexy,p,'linear');
 np=size(obj.Segments,1);
 % Special care to the last point of the crack
 L=minind==np;
 if any(L)
     minind(L)=np-1;
 end
 x1=obj.Segments(minind,2);
 y1=obj.Segments(minind,3);
 x2=obj.Segments(minind+1,2);
 y2=obj.Segments(minind+1,3);
 % sign(intvec*crackvec)
 intx1=pnts(1,1); inty1=pnts(1,2); intx2=pnts(end,1); inty2=pnts(end,2);
 signd=sign((x2-x1)*(intx2-intx1)+(y2-y1)*(inty2-inty1));
 % make the intersections points aligned with crack direction as from 1st
 % tip to 2nd tip.
 if signd==1
    pnts=flip(pnts,1);      % flip upside down
    localpnts=flip(localpnts,1);
 end

%% Check if this element is a good element (without extreme small subpolygon)
% if there is a repeated node in the pairs, meaning that the crack crosses
% two connected edges of the element, we just need calculate one triangle area
% find the repeated vertex
uniquenode=unique(pairs);
ncount=histc(pairs,uniquenode);
Lrep=ncount==2;
% the node id within the element to form a triangle with the two intersections
idrep=pairs(Lrep);
if ~isempty(idrep)
    coords=[elem.X(idrep),elem.Y(idrep)];
    l1=sqrt((coords(1)-pnts(1,1))^2+(coords(2)-pnts(1,2))^2);
    l2=sqrt((coords(1)-pnts(2,1))^2+(coords(2)-pnts(2,2))^2);
    l3=sqrt((pnts(1,1)-pnts(2,1))^2+(pnts(1,2)-pnts(2,2))^2);
    s1=(l1+l2+l3)./2;
    area=real(sqrt(s1.*(s1-l1).*(s1-l2).*(s1-l3)));
    if area/elem.Area<1e-3
        good=false;
        return;
    end
else
    % otherwise,calculate the area of the four triangles formed by the two
    % intersections and each node
    pseudoareas=zeros(1,2);
    ps1=R(1)+(1-R(2));
    ps2=R(2)+(1-R(1));
    pseudoareas(1)=ps1/ps2;
    pseudoareas(2)=ps2/ps1;
    if any(pseudoareas<2e-3)
        warning('The enriched elements will not be connected')
        good=false;
        return
    end
end
%% Select enough internal points to discretize the element 03282019
% % midpoint of the intersecption (old method was too coarse)
% midpoint=1/2*localpnts(1,:)+1/2*localpnts(2,:);
% 03282019 Try to discretize the two subpolygon finer 

%% Now move the seeds generation into a separate function: plane_partition
% scenario 1: a triangle and a pentagon--> four triangles
% if ~isempty(idrep)
%     temp=[node;localpnts];
%     tri=delaunay(temp(:,1),temp(:,2));
%     seeds=cell(1,size(tri,1));
%     for itri=1:size(tri,1)
%        id=tri(itri,:);
%        vertices=temp(id,:);
%        cent=sum(vertices)/3;    % centroid
%        seeds{itri}=[vertices;cent];% seeds for each triangle
%     end
% else    % secenario 2: two quadrilaterals
%     % subpolygons (the second and the third in the pairs must be on the same side of the crack)
%     nodesid1=pairs([1,4]);
%     nodesid2=pairs([2,3]);
%     nodes1=node(nodesid1,:);
%     nodes2=node(nodesid2,:);
%     midpoint1=sum(nodes1)/2;
%     midpoint2=sum(nodes2)/2;
%     midpointc=sum(localpnts)/2; % sum by column to give the midpoint of the crack
%     midpoint1=sum([midpoint1;midpointc])/2;
%     midpoint2=sum([midpoint2;midpointc])/2;
%     seeds1=[nodes1;localpnts;midpoint1;midpointc];
%     seeds2=[nodes2;localpnts;midpoint2;midpointc];
%     seeds={seeds1,seeds2};
% end
% seeds are the vertices for delaunay to create triangular subdomains
end

