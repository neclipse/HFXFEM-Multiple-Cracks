function BCType = CorrectBCTable_v6(Boundedge,VX,VY,BCType,fd,fd_num,BCcode)

% function BCType = CorrectBCTable(EToV,BCType,fd,BCcode);
% Purpose: Setup BCType for boundary conditions in 2D
%
%    EToV   : Element-To-Vertice table
%    VX, VY : (x,y)-coordinates of mesh vertices
%    BCType : Table with types of faces for BC's
%    fd     : handle to distance function
%    BCcode : Integer for specific boundary type
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
pxc = 0.5*(px1+px2);    % x-coordinate of edge centers
pyc = 0.5*(py1+py2);    % y-coordinate of edge centers
% calculating edge length
el=sqrt((px1(:)-px2(:)).^2+(py1(:)-py2(:)).^2);
dc = abs(fd([pxc(:) pyc(:)])); % distances to boundaries from face centers
tol=el*1E-1; % tolerance
ind= find(dc<tol);
node=Boundedge(ind,:);
fd_num = fd_num*ones(length(ind),1);
BCcode = repmat(BCcode,length(ind),1);
BCType=[BCType;fd_num node BCcode]; % Comprehensive Boundary condition table
end