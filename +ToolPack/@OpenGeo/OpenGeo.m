classdef OpenGeo < ToolPack.Geometry
    properties (Access = public)
        Id
        Linetype                                                             % 1. piecewise linear Segments, 2. function handle with two limits
        Minelength                                                           % minimum element length among elements interacted with the crack
        Maxelength
        Length                                                               % crack length calculated by callength
        Segments                                                             % n*3 matrix, each row represent one point to describe the open geometry, the larger n the better for type 1
        Description															 % description of the initial crack using either piecewise points or function handle in terms of x or y
        Boundhandles                                                         % The distance handle of the boundaries
        Phi                                                                  % The signed distance function to crack body of the obj.Nodes
        Nodespool                                                            % All the nodes within the narrow band
        Phipool                                                              % the signed disctance function to crack body of all nodes within the narrow band
        Intelements                                                          % All Enriched elements={bodyelements,rtipelements}
        Bodyelements                                                         % indices of the crack body elements
        Rtipelements                                                         % indices of the tip elements
        Blendingelems                                                        % Elements around the enriched elements meeting certain criterion, for different permeability settings.
        Nodes                                                                % All interacted within the interacted elements
        Bodynodes                                                            % All Nodes in Bodyelements
        Rtipnodes                                                            % All Nodes in Rtipelements 
        Stdnodes                                                             % A pool of nodes which are lately decided not to be enriched although they are within obj.Nodes
    end
    properties (Access = protected)
        ACTIVETIPS                                                           % describe the active Ends (notified by enrcracktip): % [1,0] the final point active and the first point inactive
    end
    properties (Dependent)
        Ends                                                                 % Two geometrical ending segements
        Omegas                                                               % tilting angle of the two Ends of the open geometry
        Tips                                                                 % Two geometrical ending points
        Rtips                                                                % Index of the Real Crack tips not lying on the geometry boundaries
    end
%     methods(Static)
%        function Id=updateid()
%           persistent counter;
%           if isempty(counter)
%               counter = 1;
%           else
%               counter=counter+1;
%           end
%           Id = counter;
%        end
%    end
    methods
        function obj = OpenGeo (id,mesh,boundhandles,nodedict,elemdict,linetype,description,ns)
            if nargin==0
                super_args={};
            else
                super_args={mesh,nodedict,elemdict};
            end
            obj = obj@ToolPack.Geometry(super_args{:});
            if nargin>3
                obj.Id=id;
                obj.Boundhandles= boundhandles;
                obj.Geotype = 1;
                obj.Linetype = linetype;
                obj.Description=description;
                % The default DIRECTION OF THE CRACK is from 2nd tip to 1st tip, 
                % i.e., the first point to the last point.
                % The normal is 90degree counterclockwisely from the crack direction.
                % Then Phi will be POSITIVE on the upper subdomain.
                % IT IS BETTER TO HAVE STRAIGHT CRACKS AT CURRENT STAGE AS
                % THE SEGMENTS DATA ARE NOT CONSIDERED IN THE INTERSECTION
                % METHOD. 11/06/2018
                if linetype==1
                    obj.Segments=description;
                    obj.discretize(ns);
                elseif linetype==2
                    obj.discretize(ns);
                end
                % description is a 3*1 cell, {1} is the function handle or segments, {2}
                % is the indicatior of independent variable, {3} is the
                % limits vector, np is the number of discrete points to be
                % generated
            end
        end
             
        function  calelen(obj)                                               % calculate the element lengths
            xlist=obj.Segments(:,2);
            ylist=obj.Segments(:,3);
            minelength=inf;
            maxelength=0;
            for ip = 1:length(xlist)
                x=xlist(ip);
                y=ylist(ip);
                [elemlist,~]=obj.Mesh.locate(x,y);                          % find connected elements to the sample point
                lengths1=obj.Elemdict(elemlist(1)).callength;               % calculate the lengths of every side of the first element
                temp1=min(lengths1);                                        % Find the shortest element length
                temp2=max(lengths1);
                if temp1<minelength                                           
                    minelength=temp1;
                end
                if temp2>maxelength
                    maxelength=temp2;
                end
            end
            obj.Maxelength=maxelength;
            obj.Minelength=minelength;                                     
        end
        
        function Ends = get.Ends(obj)
            ends1=[obj.Segments(end-1,2:3);obj.Segments(end,2:3)];               % The final Segments is defined as the first end, nodes ordered along the crack grow direction
            ends2=[obj.Segments(2,2:3);obj.Segments(1,2:3)];                     % The first Segments is defined as the last end, nodes ordered along the crack grow direction
            Ends={ends1;ends2};
        end
        
        function Tips = get.Tips(obj)
            tip1=obj.Segments(end,2:3);                                        % The last point is defined as the first tip
            tip2=obj.Segments(1,2:3);                                          % The first point is defined as the last tip
            Tips=[tip1;tip2];                                                   %  containing two tips
        end
        % Give the real crack tips not lying on the boundaries of the
        % domain
        function Rtips = get.Rtips(obj)
            % This function need to updated before being called in
            % EnrCrackTip.update. Alternatively, the ACTIVETIPS should be
            % defined to detect the active tips that are able to propagate.
            % Then we can leave the Rtips as tip not reaching the
            % boundaries.
            tip1=obj.Segments(end,2:3);                                       % The last point is defined as the first tip
            tip2=obj.Segments(1,2:3);                                         % The first point is defined as the last tip
            table=ones(length(obj.Boundhandles),2);                           % Table to store the distance of tips to the boundaries
            for ib=1: length(obj.Boundhandles)
                fb=obj.Boundhandles{ib};
                table(ib,:)=[fb(tip1),fb(tip2)];
            end
            logtable=(abs(table)>1e-9);                                     % logical table showing how close the tip is close to a boundary
            logtable=all(logtable,1);                                       % all (logtable,1) will test if the tips are not on any edge, one column one tip
            Rtips=find(logtable);                                           % index of the real tips
        end
        
        function Omegas = get.Omegas(obj)
            angle=zeros(1,2);
            CRACK = obj.Segments(:,2:3);
            disc  = CRACK(end,:)-CRACK(end-1,:);                            % Horizontal and vertical distances for final crack segment
            angle(1) = atan2(disc(2),disc(1));                              % Final Crack angle with respect to horizontal
            disc  = CRACK(1,:)-CRACK(2,:);                                  % Horizontal and vertical distances for first crack segment
            angle(2) = atan2(disc(2),disc(1));                              % First Crack angle with respect to horizontal
            Omegas=angle;                                                   % Crack angle of the geometrical crack tips
        end
        
        function [xyl,r]=callocal(obj,ict,x,y)                              % calculate the local coordinates (xl,yl) of a point (x,y), and the distance from the ict tip
            angle=obj.omega(ict);
            xct=obj.Tips{ict}(1);
            yct=obj.Tips{ict}(2);
            T=[cos(angle) sin(angle); -sin(angle) cos(angle)];              % transformation matrix
            xc=x-xct;
            yc=y-yct;
            xyl=T*[xc;yc];
            r=norm(xyl);
        end
        
        function findblending(obj,mode,varargin)
            % use a ellipse to determine elements around the injection
            % point to make the isolation zone to concentrate fluid
            % injection.
            switch mode
                % called by geo.initiate
                case 1
                    shortaxis=obj.Length;
                    injectionpoint=varargin{1};
                    if nargin==4
                        longaxis=varargin{2}*shortaxis;
                    else
                        longaxis=shortaxis;
                    end
                    x=injectionpoint(1);
                    y=injectionpoint(2);
                    X=obj.Mesh.VX;
                    Y=obj.Mesh.VY;
                    xdis=X-x;
                    ydis=Y-y;
                    dis=(xdis/shortaxis).^2+(ydis/longaxis).^2;
                    nodeid=find(dis<1);
                    elems=obj.Mesh.nodetoelem(nodeid);
                    obj.Blendingelems=elems;
                case 2
                    % called by geo.update, updated on 020920
                    if ~isempty(varargin)
                        depth=varargin{1};
                    else
                        depth=1;
                    end
                    nodeid=obj.Nodes;
                    for id=1:depth
                        elems=obj.Mesh.nodetoelem(nodeid);
                        nodeid=obj.Mesh.EToV(elems,:);
                        nodeid=unique(nodeid(:));
                    end
                    elems=setdiff(elems,obj.Intelements);
                    obj.Blendingelems=unique([obj.Blendingelems;elems]);
            end
        end
        initiate(obj,injectionpoint,varargin)                                          % to find the inititial interacted elements with the geometry
        discretize(obj,np)
        rtipelems=Toedge1(obj)                                               % Ensure the current tip lies on the edge 
        rtipelems=Toedge2(obj,varargin)                                      % Ensure the current tip lies on the edge
        update(obj)                                                         % update the geometry and the interacted elements with the updated geometry
        flag=interactwith(obj,elem)                                         % determine if the element is interacting with the crack using the calculated level set values
        length=callength(obj)                                               % determine the crack length
        h=plotme(obj,deformflag,crackflag, nodeflag,varargin);                % plot the crack segments with options
        elems=findelems( obj,plist,varargin )                               % find the interacting elements with a list of plist(x,y;...) with options 'inside', 'edge', 'in_edge'
        [good,pnts,localpnts]=intersection(obj,elem,varargin)                  % the intersection method to find the intersections of the crack with the elem obj
    end
end