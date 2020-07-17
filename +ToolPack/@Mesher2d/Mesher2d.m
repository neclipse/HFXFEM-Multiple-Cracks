% This basic class is to generate a satisfied mesh for the domain geometry
% The whole procedure is made up of three main steps:
% 1. Initially mesh the domain to give nodes and element information with
% distmesh toolbox
% 2. Validate and reorder the initial element vertices
% 3. Setup the Boundary condition table.
classdef Mesher2d < handle
   properties
   % Domain geometry

   % first meshinfo
   p         %[x,y,...] coordinates of every nodes
   t         %[ind1,ind2,ind3,...] element connections 
   % t is ordered in counterclockwise by vertices first and side nodes last 
   q         % mesh quality 
   MT        % material property type of elements
   VT        % type of mesh vertices: 1--corner nodes; 2--side nodes
   end
   properties(Dependent)
      %% second meshinfo 
      % Just to keep the compatibility with other open-source codes
      EToV      % Element to vertices, same as t.    
      VX        % x-coordinates of mesh vertices
      VY        % y-coordinates of mesh vertices
      Totelem   % total number of elements
      Totnodes  % total number of nodes
      Nface     % number of faces of an element
   end
   
   methods
       function obj=Mesher2d(p,t)
          obj.p=p;
          obj.t=t;
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
       
       function VX=get.VX(obj)
          VX=obj.p(:,1); 
       end
       
       function VY=get.VY(obj)
          VY=obj.p(:,2); 
       end
       function nodetype(obj)
          obj.VT=ones(size(obj.p,1),1);
       end      
       
% 	   function id=itemize(obj,x,y) % not good
%             % This is a method of Mesher2d to retrieve the index of a node
%             tol=1e-6;
%             id=[0,0];                                % arbitrary condition to go to while loop
%             while ~isscalar(id)                      % id is not a scalar: multiple points have been found
%             ind= find(abs(obj.VX-x)<1e-5);           % Find the index of nodes with VX==x
%             indt=logical(abs(obj.VY(ind)-y)<1e-5);   % Search within the found index for nodes with VY==y
%             id=ind(indt);                            % Convert logical value to original index
%             tol=tol/10;
%             end
%        end
       function [elemnum,nodenum,flagn] = locate(obj,x,y)                         % locate one point(x,y) inside the mesh, give the closest node and all connected elements
           X=obj.VX;
           Y=obj.VY;
           xdis=X-x;
           ydis=Y-y;
           dis=sqrt(xdis.^2+ydis.^2);                                       % dis is a list of displacements of all nodes to the given point 
           [mindis,nodenum]= min(dis);                                           % mindis is the minimun displacement, nodenum is the index of the corresponding node
           if mindis< 1e-9
               flagn=true;                                                  % determine the point almost coincides with a node
           else
               flagn=flase;
           end
           % locate connected elements with numnode
           check=logical(obj.EToV==nodenum);                                % compare the Element to vertices file to the found nodenum
           [elemnum,~]=find(check);                                         % obtain the indices of the elements connected to the nodenum 
       end
   end
    
end
