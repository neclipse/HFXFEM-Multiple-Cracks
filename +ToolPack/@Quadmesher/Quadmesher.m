classdef Quadmesher < handle
    %% a specifically designed mesher for the quarter wellbore model
    properties
        p			% nodes coordinates (x1,y1;x2,y2;...)
        t       	% connection (counterclockwisely ordered)
        VT			% nodestype defined by obj.nodestype
        radnodes	% nodes number along radial direction/width
        tannodes	% nodes number along tangential direction/height
        Nanfield  	% zeros for the mesh plot
        VX
        VY
    end
    properties(Dependent)
        %%- second meshinfo
        % Just to keep the compatibility with other open-source codes
        EToV      % Element to vertices, same as t.
        Totelem   % total number of elements
        Totnodes  % total number of nodes
        Nface     % number of faces of an element
    end
    methods
        %Constructor
        %         [ node,element,~,~ ] = mesh( re,rw,partition,tanstep,radstep1,radstep2,bias1,bias2);
        function obj=Quadmesher(node,element,varargin)
            obj.p=node;
            obj.t=element;
            obj.VX=node(:,1);
            obj.VY=node(:,2);
            obj.Nanfield=zeros(length(node),1);
            obj.reorder; %Added on 020920
            inp=inputParser;
            addOptional(inp,'radnodes',0);
            addOptional(inp,'tannodes',0);
            parse(inp,varargin{:});
            obj.radnodes=inp.Results.radnodes;
            obj.tannodes=inp.Results.tannodes;
            
        end
        function reorder(obj)
           % 020920 Make the local nodes order consistent with parent element
           for ielem=1:obj.Totelem
               nodes=obj.t(ielem,:);
               try
                    coord=obj.p(nodes,:); % all coordnates of all nodes
               catch ME
                   disp(ielem);
                   if (strcmp(ME.identifier,'MATLAB:catenate:dimensionMismatch'))
                       msg = ['Dimension mismatch occurred: First argument has ', ...
                           num2str(size(A,2)),' columns while second has ', ...
                           num2str(size(B,2)),' columns.'];
                       causeException = MException('MATLAB:myCode:dimensions',msg);
                       ME = addCause(ME,causeException);
                   end
                   rethow(ME)
               end
               
               B=sum(coord,1)/obj.Nface;
               rel_coord=coord-B;
               iso_coord=sign(rel_coord);
               % Make sure it is coutclockiwisely ordered and started from
               % the leftmost node
               n1=find(iso_coord(:,1)==-1&iso_coord(:,2)==-1);
%                n2=find(iso_coord(:,1)==1&iso_coord(:,2)==-1);
%                n3=find(iso_coord(:,1)==1&iso_coord(:,2)==1);
%                n4=find(iso_coord(:,1)==-1&iso_coord(:,2)==1);
               index1=n1:1:4;
               index2=1:(4-length(index1));
               index=[index1,index2];
               try
                    obj.t(ielem,:)=nodes(index);
               catch ME
                   if (strcmp(ME.identifier,'MATLAB:catenate:dimensionMismatch'))
                       msg = ['Dimension mismatch occurred: First argument has ', ...
                           num2str(size(A,2)),' columns while second has ', ...
                           num2str(size(B,2)),' columns.'];
                       causeException = MException('MATLAB:myCode:dimensions',msg);
                       ME = addCause(ME,causeException);
                   end
                   rethow(ME)
               end
                   
               
           end
        end
        
        function Totelem=get.Totelem(obj)
            Totelem=size(obj.t,1);
        end
        
        function Nface=get.Nface(obj)
            Nface=size(obj.t,2);
        end
        
        function Totnodes=get.Totnodes(obj)
            Totnodes=size(obj.p,1);
        end
        
        function EToV=get.EToV(obj)
            EToV=obj.t;
        end
        function nodetype(obj)
            obj.VT=ones(size(obj.p,1),1);
        end
        
        function [elemnum,nodenum,flagn] = locate(obj,x,y)                   % locate one point(x,y) inside the mesh, give the closest node and all connected elements
            X=obj.VX;
            Y=obj.VY;
            xdis=X-x;
            ydis=Y-y;
            dis=sqrt(xdis.^2+ydis.^2);                                       % dis is a list of displacements of all nodes to the given point
            [mindis,nodenum]= min(dis);                                      % mindis is the minimun displacement, nodenum is the index of the corresponding node
            if mindis< 1e-9
                flagn=true;                                                  % determine the point almost coincides with a node
            else
                flagn=false;
            end
            % locate connected elements with numnode
            % elemnum=obj.nodetoelem(nodenum)
            check=(obj.EToV==nodenum);                                % compare the Element to vertices file to the found nodenum
            [elemnum,~]=find(check);                                         % obtain the indices of the elements connected to the nodenum
        end
        function elemnum=nodetoelem(obj,nodeid)
            % give interacted elements with the given nodeid based on EToV
            nnode=numel(nodeid);
            elemnum=zeros(4*nnode,1);
            icount=1;
            for in=1:nnode
                check=logical(obj.EToV==nodeid(in));
                [elems,~]=find(check);
                icount2=icount+length(elems)-1;
                elemnum(icount:icount2)=elems;
                icount=icount2+1;
            end
            elemnum=elemnum(1:icount2);
            elemnum=unique(elemnum);
        end
        plotmesh(obj,varargin);
        elems=findelems(obj,plist,elemdict,varargin);                        % find the interacting elements with a list of plist(x,y;...) with options 'inside', 'edge', 'in_edge'
        
    end
end
