function appbound( obj )
%APPBOUND This function applies essential boundary condition to the
%linear system
%% Apply boundary condition from BCTable_line
nonodes=obj.ElemDict(1).NoNodes; %   nonodes: the number of nodes invloved in the element stiffness matrix
[~,Nline_xi]=lineshape(nonodes,1);
n=length(Nline_xi); % A indirect way to decide the number of nodes along a side of element
indx=logical(obj.BCTablel(:,2+n)==1); % Column 2+n refers to the dof of x displacement
indy=logical(obj.BCTablel(:,3+n)==1); % Column 3+n refers to the dof of y displacement
xbctable=obj.BCTablel(indx,:);        % BCTable with X-disp specified
ybctable=obj.BCTablel(indy,:);        % BCTable with Y-disp specified
locarrayX=zeros(1,n*size(xbctable,1));
locarrayY=zeros(1,n*size(ybctable,1));
% for loop is used because multilevel indexing is forbidden for vector
nx=1;ny=1;
for i=1:size(xbctable,1) % create the location array side by side, node by node
    for j=2:1+n
        locarrayX(nx)=obj.NodeDict(xbctable(i,j)).DofArray(1);
        nx=nx+1;
    end
end
for i=1:size(ybctable,1) % create the location array side by side, node by node
    for j=2:1+n
        locarrayY(ny)=obj.NodeDict(ybctable(i,j)).DofArray(2);
        ny=ny+1;
    end
end
locarray=[locarrayX,locarrayY];
%% Apply boundary condition from BCTable_node
if ~isempty(obj.BCTablen)
    xlistt=logical(obj.BCTablen(:,2+n));            % Index list of nodes with x-disp condition
    ylistt=logical(obj.BCTablen(:,3+n));            % Index list of nodes with y-disp condition
    xlist=obj.BCTablen(xlistt,1);
    ylist=obj.BCTablen(ylistt,1);
    locarraynxt=[obj.NodeDict(xlist).DofArray];     % Global x-y-dof index of specified nodes
    locarraynx=locarraynxt(1:2:end-1);              % Global x-dof index of specified nodes
    locarraynyt=[obj.NodeDict(ylist).DofArray];       % Global x-y-dof index of specified nodes
    locarrayny=locarraynyt(2:2:end);                % Global y-dof index of specified nodes
    locarray=[locarrayX,locarrayY,locarraynx,locarrayny];
end
obj.LHS(locarray,:)=0;
obj.LHS(:,locarray)=0;
obj.RHS(locarray)=0;
index = sub2ind(size(obj.LHS),locarray ,locarray );
obj.LHS(index)=1;
end
