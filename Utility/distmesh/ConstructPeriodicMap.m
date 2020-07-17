function [mapM,mapP] = ConstructPeriodicMap(EToV,VX,VY,BCType,BCcode1,BCcode2,fd)
%
% Purpose: construct periodic maps for a square.
%
% Assumptions:
%   - element nodes are ordered counterclockwise on the elements
%   - mesh is conforming at the periodic boundary
%   - distance function is a line which is angled to both boundaries
%     and is used for sorting the boundary faces in ascending/desending
%     orders
% 
% By Allan P. Engsig-Karup.

Globals2D;

VNUM = [1 2;2 3;3 1]; % face orientations

pxc = 0.5*(VX(EToV)+VX(EToV(:,[2 3 1])));
pyc = 0.5*(VY(EToV)+VY(EToV(:,[2 3 1])));
% distances to a boundary from face centers using a distance function
% the boundary defined by the distance function is used for sorting
% and subsequent matching of the faces.
dc  = abs(fd([pxc(:) pyc(:)])); 

% fetch indexes for boundaries that are to be matched
[eids1,fids1] = find(BCType == BCcode1);
[eids2,fids2] = find(BCType == BCcode2);
[idx1] = find(BCType == BCcode1);
[idx2] = find(BCType == BCcode2);
[dc1,i1] = sort(dc(idx1),'ascend');
[dc2,i2] = sort(dc(idx2),'descend');
% idx1 = idx1(i1);
% idx2 = idx2(i2);
% eids1 = eids1(i1); fids1 = fids1(i1);
eids2 = eids2(i2); fids2 = fids2(i2);
eids1 = eids1(i1); fids1 = fids1(i1);

% initialize length of new map
mapM = [];
mapP = [];

% use that on opposing boundaries the node ordering will also be opposite
% when the nodes are ordered counterclockwise on the elements
% node-ordering reversed for face 3 compared to face 1 and 2
for n = 1 : length(fids1) % go over each boundary face of BCcode type
    if fids1(n)==3
        mapM = [mapM (eids1(n)-1)*Nfaces*Nfp+[ (fids1(n)-1)*Nfp + [Nfp:-1:1] ]];
    else
        mapM = [mapM (eids1(n)-1)*Nfaces*Nfp+[ (fids1(n)-1)*Nfp + [1:Nfp] ]];
    end
    % note node orderings reversed on opposing boundaries...
    if fids2(n)==3
        mapP = [mapP (eids2(n)-1)*Nfaces*Nfp+[ (fids2(n)-1)*Nfp + [1:Nfp] ]];
    else
        mapP = [mapP (eids2(n)-1)*Nfaces*Nfp+[ (fids2(n)-1)*Nfp + [Nfp:-1:1] ]];
    end
end

%% Visual sanity check
% figure
% for i = 1 : length(mapM)
%     clf
%     triplot(EToV,VX,VY,'k')
%     hold on
%     plot(x,y,'k.','markersize',20)
%     plot(Fx(mapM(i)),Fy(mapM(i)),'ro')
%     plot(Fx(mapP(i)),Fy(mapP(i)),'bo')
%     drawnow
%     pause
% end
