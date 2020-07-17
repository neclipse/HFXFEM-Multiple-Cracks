function initialRHS( obj )
%CRTRHS_V3 Create the right hand side of the linear system
%   This function include integration and assembly
%   This new version is compatible with new constrBCTable.m for two phase
%   problem.04242018
agl=obj.ElemDict(1).GaussPntDictM(1).GBINP;
inistress=agl.inistress;
INI=zeros(obj.Drow,1);
nonodes=obj.ElemDict(1).NoNodes;                    % nonodes: the number of nodes invloved in the element stiffness matrix
[~,Nline_xi]=lineshape(nonodes,1);
n=length(Nline_xi);                                 % An indirect way to decide the number of nodes along a side of element
if nonodes==3 || nonodes==6                         % For triangular elements
    [xx,ww]=GaussQuad(1,1);                         % calculating gauss points along a line from 0 to 1
elseif nonodes==4 || nonodes==8                     % For rectangular elements
    [xx,ww]=GaussQuad(1);                           % calculating gauss points along a line from -1 to 1
end
%% -- Calculate the initial load vector
if ~isempty(obj.BCTableim)  % may not need to calculate obj.INI if external bounaries are all fixed, 01/23/19
    table=obj.BCTableim;
    for i=1:size(table,1)
        fe=zeros(2*n,1);
        locarray=zeros(2*n,1);
        X=[obj.NodeDict(table(i,1:n)).X];             % X is a row vector
        % X is the x coordinates list of the nodes along the boundary side
        Y=[obj.NodeDict(table(i,1:n)).Y];
        % loop over gaussian points
        for j=1:length(xx)
            xi=xx(j);
            [Nline,Nline_xi]=lineshape(nonodes,xi);
            x_xi=Nline_xi*X';
            y_xi=Nline_xi*Y';
            normal=[y_xi,0,-x_xi,0;0,-x_xi,y_xi,0];    % THIS FORMULA IS DEBATABLE for Triangle Elements
            fet=Nline'*normal*inistress;               % THIS FORMULA IS DEBATABLE for Triangle Elements
            fe=fe+ww(j)*fet;
        end
        %%- Assembly into global force vector
        for k=1:n
            %IMPORTANT ERROR: TABLE(I,K+1) WAS CODED AS TABLE(I,K)--050918
            locarray(2*k-1:2*k)=obj.NodeDict(table(i,k)).DofArray(1:2);
            % using doflist in the node object to construct a locarray for
            % assembly
        end
        INI(locarray)=INI(locarray)+fe;
    end
end
obj.INI=INI;
end                              


