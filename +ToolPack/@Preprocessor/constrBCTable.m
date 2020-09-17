function constrBCTable(obj)
%% Construct the BCTable_line first, this is the basis of boundary condition table
%RHS are buided directly based on BCTable_line
%BCTable_disp is a combination of BCTable_line and BCTable_node. It is builded for appbound.m
%Special nodes(which have been already constrained by boundary loads)
%should be avoided in the BCTable_disp.
tol=1E-5;
maxiter=3;
iiter=0;
nln=size(obj.BoundEdge,2); % number of nodes on a side of an element
while size(obj.BCTable_line,1)~=size(obj.BoundEdge,1)
    obj.BCTable_line=[];   % set a void BCTable_line
    for i=1:obj.Nobound
        fd=obj.Disthandle{i};
        fd_num=i;
        BCcode=obj.BCmat_line(i,:);
        obj.BCTable_line=CBCTable(obj.BoundEdge,obj.Mesher.VX,obj.Mesher.VY,obj.BCTable_line,fd,fd_num,BCcode,tol);
        % fd: handle of the distance function of a boundary 
        % BCcode: integer vertice representing the boundary condition
        % of all Dofs. For example, [1,1,1] means constant displacement boundary
        % condition applies to the first two Dofs and constant load
        % condition applies to the third Dof.
        % tol: the numerical tolerance of distance to locate the boundeges
    end
    if size(obj.BCTable_line,1)>size(obj.BoundEdge,1)
        tol=tol*10;
        elseif size(obj.BCTable_line,1)<size(obj.BoundEdge,1)
        tol=tol/10;
    end
    iiter=iiter+1;
    if iiter>maxiter
       error('BCTable_line is built with error.') 
    end
end
%% construct the imbalanced external boundary for linsyscreater.initialRHS
% Neumann boundary condition
imbalance=[];
if ~isempty(obj.Imbalanced)
    for im=1:length(obj.Imbalanced)
        idm=obj.Imbalanced(im);
        [irm,~]=find(obj.BCTable_line(:,1)==idm);
        imbalance=[imbalance;obj.BCTable_line(irm,2:3)];
    end
end
obj.BCTable_imbalance=imbalance;
%% Construct the BCTable_node accordingly
if ~isempty(obj.BCTable_node)
    for i=1:size(obj.BCTable_node,1)
        % The specified nodes are given coordinates or index. Index for 
        % those with coordinates is preassigned as zero, and should be
        % retrieved through the "itemize" function
        if obj.BCTable_node(i,2)==0
           [~,nodeid]=obj.Mesher.locate(obj.BCTable_node(i,3),obj.BCTable_node(i,4));
           obj.BCTable_node(i,2)=nodeid;
        end
    end
    idnode=obj.BCTable_node(:,2);               % Index of all special(contradictory) nodes
    % Format of BCTable_node:
    % 1st column: The type of the dof, 1- xdisp, 2-ydisp, 3-pressure, 4,
    % xdisp_enr, 5-ydisp_enr
    % 2nd column: The global index of the node
    % 3rd and 4th: The coordinates of the node
    % 5th column: The prescribed displacement or pressure value
else
    idnode=[];
end
%% Construct BCTable_Dirichlet
[ir1,ic1]= find(obj.BCmat_line==1);
obj.BCTable_disp=cell(size(ir1,1),1);
for idisp=1:length(ir1)
   indbd=logical(obj.BCTable_line(:,1)==ir1(idisp)); % Find all edges on this boundary
   bctemp=obj.BCTable_line(indbd,:);             
   bctemp=bctemp(:,2:nln+1);
   bctemp=unique(bctemp);
   bctemp=bctemp';
   if ~isempty(idnode)                              % If some nodes are already prescribed in the bctable-node, exlude them here
       for idn=1:length(idnode)                     % Delete speical nodes one by one
          bctemp=bctemp(bctemp~=idnode(idn)); 
       end
   end
   idof = repmat(ic1(idisp),length(bctemp),1);       % boundary condition, 1-x_disp; 2-y-disp, 3-p
   dof=repmat(obj.Predof(ir1(idisp),ic1(idisp)),length(bctemp),1);   % predefined dof value;
   bcdiri=[bctemp;idof';dof'];
   obj.BCTable_disp{idisp}=bcdiri;
end

%% Construct load boundary condition table BCTable_Neumann
[ir2,~]= find(obj.BCmat_line==2);                           % Stress boundary load
[ir3,~]= find(obj.BCmat_line==3 | obj.BCmat_line== 4);      % Traction boundary load
[ir4,~]= find(obj.BCmat_line==5 | obj.BCmat_line== 6);      % flow boundary load
obj.BCTable_stress=cell(size(ir2,1),1);
obj.BCTable_traction=cell(size(ir3,1),1);
obj.BCTable_flow=cell(size(ir4,1),1);
    %-- Build BCTable_stress for boundaries with uniform spreading load in
    % cartesian coordinate system
    for icart=1:length(ir2)                                 
        indbd=logical(obj.BCTable_line(:,1)==ir2(icart)); % Find all edges on this boundary
        bctemp=obj.BCTable_line(indbd,:);
        bctemp=bctemp(:,2:nln+1);
        curloads=obj.CurLoads{ir2(icart),1};
        Curloads=repmat(curloads',size(bctemp,1),1);
        bcstress=[bctemp,Curloads];
        obj.BCTable_stress{icart}=bcstress;
    end

    %-- Build BCTable_traction for boundaries with varing spreading load in
    % cartesian coordinate system. [tn,ts]-->[tx,ty]
    % However, the method to apply traction as boundary condition in crtRHS is not validated yet.
    for icylind=1:length(ir3)
        indbd=logical(obj.BCTable_line(:,1)==ir3(icylind));   % Find all edges on this boundary
        bctemp=obj.BCTable_line(indbd,:);
        bctemp=bctemp(:,2:nln+1);
        curloads=obj.CurLoads{ir3(icylind),1};
        if obj.BCmat_line(ir3(icylind),4)==3
            if nln==2
                x1=obj.Mesher.p(bctemp(:,1),1);                % x-coord of the first node
                x2=obj.Mesher.p(bctemp(:,2),1);                % x-coord of the second node
                X=(x1+x2)/2;
                y1=obj.Mesher.p(bctemp(:,1),2);                % y-coord of the first node
                y2=obj.Mesher.p(bctemp(:,2),2);                % y-coord of the second node
                Y=(y1+y2)/2;
            elseif nln==3
                X=obj.Mesher.p(bctemp(:,2),1);
                Y=obj.Mesher.p(bctemp(:,2),2);
            end
            % Transform stress form cylindrical coordinate system to cartesian coordinate system
            Curloads=loadr2x(curloads,X,Y);
        elseif obj.BCmat_line(ir3(icylind),4)==4
            Curloads=repmat(curloads,size(bctemp,1),1);
        end
        bctraction=[bctemp,Curloads];
        obj.BCTable_traction{icylind}=bctraction;
    end
    
    % -- Build the volumetric flow condition using flow vector [qx,qy]-->qn
    for icylind=1:length(ir4)                                 
        indbd=logical(obj.BCTable_line(:,1)==ir4(icylind));   % Find all edges on this boundary
        bctemp=obj.BCTable_line(indbd,:);
        bctemp=bctemp(:,2:nln+1);
        curloads=obj.CurLoads{ir4(icylind),2};
        if obj.BCmat_line(ir4(icylind),3)==5
            Curloads=repmat(curloads,size(bctemp,1),1);
        elseif obj.BCmat_line(ir4(icylind),3)==6
            if nln==2
                x1=obj.Mesher.p(bctemp(:,1),1);                % x-coord of the first node
                x2=obj.Mesher.p(bctemp(:,2),1);                % x-coord of the second node
                X=(x1+x2)/2;
                y1=obj.Mesher.p(bctemp(:,1),2);                % y-coord of the first node
                y2=obj.Mesher.p(bctemp(:,2),2);                % y-coord of the second node
                Y=(y1+y2)/2;
            elseif nln==3
                X=obj.Mesher.p(bctemp(:,2),1);
                Y=obj.Mesher.p(bctemp(:,2),2);
            end
%       Transform x-y form flow vector to normal flow and shear flow
            Curloads=loadr2r(curloads,X,Y);
        end
        bcflow=[bctemp,Curloads(:,1)];                     % only store normal component of the flow
        obj.BCTable_flow{icylind}=bcflow;
    end
end




















