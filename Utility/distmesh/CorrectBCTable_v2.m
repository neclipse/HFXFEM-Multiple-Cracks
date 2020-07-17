function BCType = CorrectBCTable_v2(EToV,VX,VY,BCType,fd,BCcode)

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
% VNUM = [1 2;2 3;3 1]; % face orientations
Globals2D;
pxc = 0.5*(VX(EToV(:,[1 2 3]))+VX(EToV(:,[2 3 1])));    % x-coordinate of face centers
pyc = 0.5*(VY(EToV(:,[1 2 3]))+VY(EToV(:,[2 3 1])));    % y-coordinate of face centers
pxc1 = VX(EToV(:,[1 2 3]));    % x-coordinate of first node of a side
pxc2 = VX(EToV(:,[2 3 1]));    % x-coordinate of second node of the side
pyc1 = VY(EToV(:,[1 2 3]));    % y-coordinate of first node of a side
pyc2 = VY(EToV(:,[2 3 1]));    % y-coordinate of second node of the side
% calculating edge length
plc=sqrt((pxc1(:)-pxc2(:)).^2+(pyc1(:)-pyc2(:)).^2);
dc  = abs(fd([pxc(:) pyc(:)])); % distances to boundaries from face centers
tol=plc*1e-1; % tolerance
idx = dc<tol;
BCType(idx) = BCcode;
end