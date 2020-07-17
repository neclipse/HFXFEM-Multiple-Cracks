function BCType = CorrectBCTable_v5(EToV,VX,VY,BCType,fd,fd_num,BCcode)

% function BCType = CorrectBCTable(EToV,BCType,fd,BCcode);
% Purpose: Setup BCType for boundary conditions in 2D
%
%    EToV   : Element-To-Vertice table
%    VX, VY : (x,y)-coordinates of mesh vertices
%    BCType : Table with types of faces for BC's
%    fd     : handle to distance function
%    BCcode : Integer for specific boundary type
%
% By Allan P. Engsig-Karup 
%
% Globals2D;
sdorder = [1 2;2 3;3 1]; % face orientations
pxc = 0.5*(VX(EToV(:,[1 2 3]))+VX(EToV(:,[2 3 1])));    % x-coordinate of face centers
pyc = 0.5*(VY(EToV(:,[1 2 3]))+VY(EToV(:,[2 3 1])));    % y-coordinate of face centers
px1 = VX(EToV(:,[1 2 3]));    % x-coordinate of first node of a side
px2 = VX(EToV(:,[2 3 1]));    % x-coordinate of second node of the side
py1 = VY(EToV(:,[1 2 3]));    % y-coordinate of first node of a side
py2 = VY(EToV(:,[2 3 1]));    % y-coordinate of second node of the side
% calculating edge length
el=sqrt((px1(:)-px2(:)).^2+(py1(:)-py2(:)).^2);
dc  = abs(fd([pxc(:) pyc(:)])); % distances to boundaries from face centers
tol=el*1e-1; % probable maximum distance from a face center to the bounday
% IMPORTANT BUG: THE TOLERANCE IS NOT CALCULATED ACCURATELY, THUS THE
% SELECTED EDGES MAY NOT BE REAL BOUNDARY EDGES
ind= find(dc<tol);
[I, J]=ind2sub(size(EToV),ind); % convert linear index to sub index to determin local face order
J2=sdorder(J,2);                % determin the local second node order
ind2=sub2ind(size(EToV),I,J2);  % return the linear index of the second node
node1=EToV(ind);
node2=EToV(ind2);
fd_num = fd_num*ones(length(I),1);
BCcode = repmat(BCcode,length(I),1);
BCType=[BCType;fd_num node1 node2 BCcode]; % Comprehensive Boundary condition table
end