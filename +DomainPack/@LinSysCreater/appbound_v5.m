function appbound_v5( obj )
%APPBOUND_v5 simplified verison of appbound_v4.
% This function applies essential boundary condition to the
%linear system using the direct method, but only handle zero condition.
% This will also work with speicial standnodes described by the
% enrichitems. 03032019

%% Loop through all dirichlet boundary condition and assemble the complete locarray of known dofs Ud
Allknownsloc=zeros(obj.Drow,1);                                            
iloc=0;
iloc0=1;
% from BCTable_disp                               
if ~isempty(obj.BCTabled)
    %     Pseudobound=[];                           
    for id=1:length(obj.BCTabled)                  
        disptable=obj.BCTabled{id};
        nodelist=disptable(1,:);
        idof=disptable(2,1);                        
        locarray=[obj.NodeDict(nodelist).DofArray]; 
        locarray=locarray(idof:3:end);
        iloc=iloc+length(locarray);
        Allknownsloc(iloc0:iloc)=locarray;
        iloc0=iloc+1;
    end
end
if ~isempty(obj.BCTablen)
    for inode=1:size(obj.BCTablen,1)
        idof=obj.BCTablen(inode,1);
        nodeid=obj.BCTablen(inode,2);
        iloc=iloc+1;
        Allknownsloc(iloc)=obj.NodeDict(nodeid).DofArray(idof);
    end
    iloc0=iloc+1;
end
if ~isempty(obj.BCTableEn)
    locarray=obj.BCTableEn(:,1);
%     knowns=obj.BCTableEn(:,2);
    iloc=iloc+length(locarray);
    Allknownsloc(iloc0:iloc)=locarray;
end
Allknownsloc=Allknownsloc(1:iloc);
[dloc,~,~]=unique(Allknownsloc);      
if length(dloc)~=iloc
%    warning('Contraditory Dirichlet boundary is encountered');
end
obj.LHS(dloc,:)=0;
obj.LHS(:,dloc)=0;
index = sub2ind(size(obj.LHS),dloc,dloc);
obj.LHS(index)=1;
obj.RHS(dloc)=0;
end


