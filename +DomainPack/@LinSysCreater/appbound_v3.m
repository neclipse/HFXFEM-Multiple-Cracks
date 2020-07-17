function appbound_v3( obj )
%APPBOUND_v3 This function applies essential boundary condition to the
%linear system. '0' dp and du.
% This version only allows incremental zero displacement boundary
% condition.
%% Apply essential boundary conditon from BCTable_disp                               
if ~isempty(obj.BCTabled)
%     Pseudobound=[];                                 % store these nodes for NewRapItr.converger
    for id=1:length(obj.BCTabled)                 % Go through all displacement conditions and pore pressure condition
       disptable=obj.BCTabled{id};
       nodelist=disptable(1,:);
       idof=disptable(2,1);                         % The idof is stored on the second row of the table
       locarray=[obj.NodeDict(nodelist).DofArray];  % Retrieve the boundary condition, 1-x_disp; 2-y-disp, 3-p
       locarray=locarray(idof:3:end);               
%        %- Approximation method to assign boundary condition
%        index = sub2ind(size(obj.LHS),locarray, locarray);
%        obj.LHS(index)=1E30;
%        obj.RHS(locarray)= 1E30*disptable(3,:);
       %- O-1 Method to assign boundary condtion
        obj.LHS(locarray,:)=0;
        obj.LHS(:,locarray)=0;
        obj.RHS(locarray)=0;
        index = sub2ind(size(obj.LHS),locarray ,locarray);
        obj.LHS(index)=1;
    end
end
%% Apply essential boundary conditon from BCTable_node (sparse nodes)
if ~isempty(obj.BCTablen)
    for inode=1:size(obj.BCTablen,1)
       idof=obj.BCTablen(inode,1);
       nodelist=obj.BCTablen(inode,2);
       locarrayn=obj.NodeDict(nodelist).DofArray(idof);
       index = sub2ind(size(obj.LHS),locarrayn ,locarrayn);
       obj.LHS(index)=1E30;
       obj.RHS(locarrayn)= 1E30*obj.BCTablen(inode,5);
    end
end

end

