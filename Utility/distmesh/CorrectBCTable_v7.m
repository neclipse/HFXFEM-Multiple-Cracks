function BCType = CorrectBCTable_v7(Boundedge,VX,VY,BCType,fd,fd_num,BCcode,tol)

% function BCType = CorrectBCTable(EToV,BCType,fd,BCcode);
% Purpose: Setup BCType for boundary conditions in 2D
%
%    EToV   : Element-To-Vertice table
%    VX, VY : (x,y)-coordinates of mesh vertices
%    BCType : Table with types of faces for BC's
%    fd     : handle to distance function
%    BCcode : Integer for specific boundary type
%    tol    : distance tolerance
%
% By Allan P. Engsig-Karup and Chang Huang
%
% Globals2D;
% sdorder = [1 2;2 3;3 1]; % face orientations
% assign bundary label to each boundary edge
px1=VX(Boundedge(:,1));    % x-coordinate of first node of a boundary edge
px2=VX(Boundedge(:,2));    % x-coordinate of second node of the boundary edge
py1=VY(Boundedge(:,1));    % y-coordinate of first node of a boundary edge
py2=VY(Boundedge(:,2));    % y-coordinate of second node of the boundary edge

dc1 = abs(fd([px1(:) py1(:)])); % distances to boundaries from first node
dc2 = abs(fd([px2(:) py2(:)])); % distances to boundaries from second node
% find edges whose both nodes are on the boundary
ind= find(dc1<tol);
indt=logical(dc2(ind)<tol);
ind=ind(indt);
nodes=Boundedge(ind,:);
fd_num = fd_num*ones(length(ind),1);
BCcode = repmat(BCcode,length(ind),1);
BCType=[BCType;fd_num nodes BCcode]; % Comprehensive Boundary condition table

end