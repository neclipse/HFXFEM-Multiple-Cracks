function appbound_v4( obj )
%APPBOUND_v4 This function applies essential boundary condition to the
%linear system using the direct method, but attempts to handle both
% zero and non-zero condition. This version will also work with the
% speicial standnodes described by the enrichitems. 03032019

%% Loop through all dirichlet boundary condition and assemble the complete locarray of known dofs Ud
Allknownsloc=zeros(obj.Drow,1);                     % locarray of all dofs described by dirichlet
Allknowns=zeros(obj.Drow,1);                        % all known dofs described by dirichlet
iloc=0;
iloc0=1;
% from BCTable_disp                               
if ~isempty(obj.BCTabled)
    %     Pseudobound=[];                            % store these nodes for NewRapItr.converger
    for id=1:length(obj.BCTabled)                    % Go through all displacement conditions and pore pressure condition
        disptable=obj.BCTabled{id};
        nodelist=disptable(1,:);
        idof=disptable(2,1);                         % The idof is stored on the second row of the table
        locarray=[obj.NodeDict(nodelist).DofArray];  % Retrieve the boundary condition, 1-x_disp; 2-y-disp, 3-p
        locarray=locarray(idof:3:end);
        iloc=iloc+length(locarray);
        Allknownsloc(iloc0:iloc)=locarray;
        Allknowns(iloc0:iloc)=disptable(3,:);
        iloc0=iloc+1;
    end
end
% from sparse nodes in obj.BCTablen
if ~isempty(obj.BCTablen)
    for inode=1:size(obj.BCTablen,1)
        idof=obj.BCTablen(inode,1);
        nodeid=obj.BCTablen(inode,2);
        iloc=iloc+1;
        Allknownsloc(iloc)=obj.NodeDict(nodeid).DofArray(idof);
        Allknowns(iloc)=obj.BCTablen(inode,5);
    end
    iloc0=iloc+1;
end
% from sparse nodes in obj.BCTableEn
if ~isempty(obj.BCTableEn)
    locarray=obj.BCTableEn(:,1);
    knowns=obj.BCTableEn(:,2);
    iloc=iloc+length(locarray);
    Allknownsloc(iloc0:iloc)=locarray;
    Allknowns(iloc0:iloc)=knowns;
end
% Delete the zeros in the two arrays
Allknowns=Allknowns(1:iloc);
Allknownsloc=Allknownsloc(1:iloc);
% Sort the locarray and Test if there are repeated dofs
[dloc,ia,~]=unique(Allknownsloc);      % the result is by default sorted
if length(dloc)~=iloc
%    warning('Contraditory Dirichlet boundary is encountered');
end
dofd=Allknowns(ia);% sorted dirichlet dofs, can be non-zero values
%% Apply essential boundary condition using the direct method by taking advantage of MATLAB submatrix
% if any(dofd)
    allloc=1:obj.Drow;                  % All dofs in the sorted order
    uloc=setdiff(allloc,dloc);          % The unknown dofs in the sorted order
    % subtract kud*dofd if dofd contains non-zeros
    kud=obj.LHSnew(uloc,dloc);
    obj.RHS(uloc)=obj.RHS(uloc)-kud*dofd;
% end
% for LHS, use "0-1" method
obj.LHSnew(dloc,:)=0;
obj.LHSnew(:,dloc)=0;
index = sub2ind(size(obj.LHSnew),dloc,dloc);
obj.LHSnew(index)=1;
obj.RHS(dloc)=dofd;
end


