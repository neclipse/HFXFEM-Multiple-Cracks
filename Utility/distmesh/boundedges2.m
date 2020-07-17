function e=boundedges2(p,t)
%BOUNDEDGES Find boundary edges from pure triangular or pure quadrilateral mesh
%01092020
%Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.
% Modified by Chang Huang, 2018-2020

% Form all edges, non-duplicates are boundary edges
if size(t,2)==3|| size(t,2)==6
    edges=[t(:,[1,2]);
           t(:,[1,3]);
           t(:,[2,3])];
    node3=[t(:,3);t(:,2);t(:,1)];
elseif size(t,2)==4||size(t,2)==8
    % the order of edges and orientation of nodes do not matter for now
    edges=[t(:,[1,2]);
           t(:,[2,3]);
           t(:,[3,4]);
           t(:,[4,1])];
    % node3 is the third node that can form a triangle with the
    % corresponding edge, i.e., [1,2]--> 3 or 4. Confirmed on 01092020.
    node3=[t(:,3);t(:,4);t(:,1);t(:,2)];
end
% sort the two nodes of each row in order to better find unique edges
edges=sort(edges,2); 
[~,ix,jx]=unique(edges,'rows');
% On 010920, change histc to histcounts by using special Integer method to
% count the occurences of each integer within jx. This method should be
% more compatible with newer version of Matlab
vec=histcounts(jx,'BinMethod','Integers'); 
% old method
% vec=histc(jx,1:max(jx));
qx=find(vec==1);
% Find the uniqe edge
e=edges(ix(qx),:);
% Find the corresponding third node to make a triangle
node3=node3(ix(qx));

% Force the node Orientation to be counter-clockwise
% vector 3->1
v1=p(e(:,1),:)-p(node3,:);
% vector 3->2
v2=p(e(:,2),:)-p(node3,:);
% The cross product of two vectors originating from node 3, positive when
% 1->2 is counter-clockwisely ordered, vice versa. 010920
ix=find(v1(:,1).*v2(:,2)-v1(:,2).*v2(:,1)<0);
e(ix,[1,2])=e(ix,[2,1]);
