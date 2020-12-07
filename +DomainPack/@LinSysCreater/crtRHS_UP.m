function crtRHS_UP( obj,newmark )
%CRTRHS_V3 Create the right hand side of the linear system
%   This function include integration and assembly
%   This new version is compatible with new constrBCTable.m for two phase
%   problem.04242018
FFG=zeros(obj.Drow,1);
nonodes=obj.ElemDict(1).NoNodes;                    % nonodes: the number of nodes invloved in the element stiffness matrix
[~,Nline_xi]=lineshape(nonodes,1);
n=length(Nline_xi);                                 % An indirect way to decide the number of nodes along a side of element
% if nonodes==3 || nonodes==6                         % For triangular elements
%     [xx,ww]=GaussQuad(1,1);                         % calculating gauss
%     points along a line from 0 to 1 (wrong)
% elseif nonodes==4 || nonodes==8                     % For rectangular elements
%(10252018-https://scicomp.stackexchange.com/questions/27441/line-integral-along-the-edge-of-an-isoparametrically-mapped-triangle)
[xx,ww]=GaussQuad(1);                               % calculating gauss points along a line from -1 to 1 regardless the element type
% end
a1=newmark.a1;                                      % Newmark parameters
%% -- Calculate the current load vector using BCTable_stress
stable=obj.BCTables;                                 % BCTable_load: node indices and side-wise stress tensor
for istress=1:size(stable,1)
    table=stable{istress};
    for i=1:size(table,1)
        fe=zeros(2*n,1);
        locarray=zeros(2*n,1);
        X=[obj.NodeDict(table(i,1:n)).X];             % X is a row vector
        % X is the x coordinates list of the nodes along the boundary side
        Y=[obj.NodeDict(table(i,1:n)).Y];
        loads=table(i,(n+1):end);
        % loop over gaussian points
        for j=1:length(xx)
            xi=xx(j);
            [Nline,Nline_xi]=lineshape(nonodes,xi);
            x_xi=Nline_xi*X';
            y_xi=Nline_xi*Y';
            normal=[y_xi,0,-x_xi,0;0,-x_xi,y_xi,0];    % THIS FORMULA IS DEBATABLE for Triangle Elements
            fet=-a1*Nline'*normal*loads';                   % THIS FORMULA IS DEBATABLE for Triangle Elements
            fe=fe+ww(j)*fet;
        end
        %%- Assembly into global force vector
        for k=1:n
            locarray(2*k-1:2*k)=obj.NodeDict(table(i,k)).DofArray(1:2);
            % using doflist in the node object to construct a locarray for
            % assembly
        end
        FFG(locarray)=FFG(locarray)+fe;
    end
end
%% -- Calculate the current load vector using BCTable_traction
ttable=obj.BCTablet;                                 % BCTable_traction: node indices and side-wise traction                              % BCTable_load: node indices and side-wise stress tensor
for itraction=1:size(ttable,1)
    table=ttable{itraction};
    for i=1:size(table,1)
        fe=zeros(2*n,1);
        locarray=zeros(2*n,1);
        X=[obj.NodeDict(table(i,1:n)).X];             % X is a row vector
        % X is the x coordinates list of the nodes along the boundary side
        Y=[obj.NodeDict(table(i,1:n)).Y];
        loads=table(i,(n+1):end);                     % loads=[tx,ty];
        % loop over gaussian points
        for j=1:length(xx)
            xi=xx(j);
            [Nline,Nline_xi]=lineshape(nonodes,xi);
            x_xi=Nline_xi*X';
            y_xi=Nline_xi*Y';
            ds=sqrt(x_xi^2+y_xi^2);
            fett=loads'*ds;
            fet=-a1*Nline'*fett;                       % THIS FORMULA IS DEBATABLE for Triangle Elements
            fe=fe+ww(j)*fet;
        end
        %%- Assembly into global force vector
        for k=1:n
            locarray(2*k-1:2*k)=obj.NodeDict(table(i,k)).DofArray(1:2);
            % using doflist in the node object to construct a locarray for
            % assembly
        end
        FFG(locarray)=FFG(locarray)+fe;
    end
end

%% -- Calculate the current load vector using BCTable_flow
ttable=obj.BCTableq;                                 % BCTable_traction: node indices and side-wise traction
for itraction=1:size(ttable,1)
    table=ttable{itraction};
    for i=1:size(table,1)
        fe=zeros(n,1);
        locarray=zeros(n,1);
        X=[obj.NodeDict(table(i,1:n)).X];             % X is a row vector
        % X is the x coordinates list of the nodes along the boundary side
        Y=[obj.NodeDict(table(i,1:n)).Y];
        loads=table(i,(n+1):end);
        % loop over gaussian points
        for j=1:length(xx)
            xi=xx(j);
            [Nline,Nline_xi]=lineshape(nonodes,xi);
            x_xi=Nline_xi*X';
            y_xi=Nline_xi*Y';
            ds=sqrt(x_xi^2+y_xi^2);
            fett=loads'*ds;
            Nlinep=[Nline(1,1),Nline(1,3)];            % Nlinep is simply 1D, [N1,N2] along the line element
            fet=Nlinep'*fett;
            fe=fe+ww(j)*fet;
        end
        %%- Assembly into global force vector
        for k=1:n
            locarray(k)=obj.NodeDict(table(i,k)).DofArray(3);
            % using doflist in the node object to construct a locarray for
            % assembly
        end
        % 04232019, CHANGE SIGN CONVENTION TO OUTWARD POSITIVE.
        FFG(locarray)=FFG(locarray)-fe;
    end
end

%% -- Calculate the current flow load vector using BCTable_flow_enr, 03172019
% on 0429, calculate both BCTable_flow_std and BCTable_flow_enr only for
% the enriched element.
if ~isempty(obj.BCTableqEn)
    for iqenr=1:length(obj.BCTableqEn)
        table=obj.BCTableqEn{iqenr};
        for iqn=1:size(table,1)
            id=table(iqn,1);
            encrack=obj.EnrichItems(id);
            q=table(iqn,2);              % averaged normal flux rate per m, in 2d the unit is m^2/s
%             nodes=table(iqn,[3,4]);      % local index of the nodes along the application edge within the enriched element.
%           no longer use nodes in table to determine the closest nodes
%           from the injection point. Instead, it is determined inside
%           cal_qextenr
            if size(table,2)>2
                injectionpoint=table(iqn,[3,4]);
                [fe,locarray]=encrack.cal_qextenr(q,injectionpoint);    % return flow volumes to the specified dofs
            else
                [fe,locarray]=encrack.cal_qextenr(q);    % default injection point is the edge crack mouth
            end
            % 04232019, CHANGE SIGN CONVENTION TO OUTWARD POSITIVE.
            FFG(locarray)=FFG(locarray)-fe;
        end
    end
end
obj.RHS=FFG-(-a1*obj.INI);                      % RHS is the change of load vector comparing to the initial stress state
end

