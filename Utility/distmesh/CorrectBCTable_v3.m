function BCType = CorrectBCTable_v3(EToV,VX,VY,BCType,fd,BCcode)

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
      
pxc1 = VX(EToV(:,[1 2 3]));    % x-coordinate of first node of a side
pxc2 = VX(EToV(:,[2 3 1]));    % x-coordinate of second node of the side
pyc1 = VY(EToV(:,[1 2 3]));    % y-coordinate of first node of a side
pyc2 = VY(EToV(:,[2 3 1]));    % y-coordinate of second node of the side
% calculating edge length
plc=sqrt((pxc1(:)-pxc2(:)).^2+(pyc1(:)-pyc2(:)).^2);
tol=plc*1e-1;
d1 = abs(fd([pxc1(:) pyc1(:)])); % distances to boundary from the first node
d2 = abs(fd([pxc2(:) pyc2(:)])); % distances to boundary from the second node
idx1=logical(d1<tol);           %  test if the first node is on the boundary
idx2=logical(d2<tol);
idx=idx1+idx2;
idx =logical(idx==2);% if the two nodes are both on the boundary, this is a boundary edge
% IMPORTANT BUG: IF THREE NODES ARE ALL ON A BOUNDARY(CURVE), THE ALGORITHM
% TAKE ALL THREE EDGES AS BOUNDARY EDGES. WRONG
                
BCType(idx) = BCcode;
end