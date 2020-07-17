classdef Preprocessor < handle
   properties
       Mesher           % an object of Mesher
       Nobound          % number of boundaries
       Disthandle       % the distance handle of the boundaries
       Psdboundind      % index the dirichlet boundary
       Imbalanced       % index of Neumann boundary
       Psdnodes          % nodes to be ignored on the pseudo-boundary
       CurLoads         % Current loads applied to the boundaries
       BCmat_line       % boundary condition matrix, its size is nobound* nodofs
       Predof
       BoundEdge        % indices of corner nodes on a side along the boundary
       BCTable_line
       BCTable_imbalance% mechanically loaded boundary (Linsys.initialRHS)
       BCTable_node     % boundary condition table, its size is nonode*4
       QTable_node      % line source/link in arbitrary location,[[localnode1,localnode2],{points},q]
       BCTable_disp     % boundary condition table for predescribed displacement
       BCTable_stress   % boundary condition table for predescribed stress tensor, TENSILE POSITIVE
       BCTable_traction % boundary condition table for predescribed traction vector
       BCTable_flow     % boundary condition table for predescribed flow vector, OUTWARD POSITIVE
       % First coloumn is the index of boundary, next two coloumns are 
       % indices of nodes on a side along the boundary, last several
       % coloumns are the boundary condition type of each kind of dof
    
   end
   
   methods
       
       function obj=Preprocessor(nobound,disthandle,BCmat_line,varargin)
           obj.Nobound=nobound;
           obj.Disthandle=disthandle;
           obj.BCmat_line=BCmat_line;
           if nargin>3
           obj.BCTable_node=varargin{1};
           obj.Predof=varargin{2};
           obj.CurLoads=varargin{3};
           obj.Psdboundind=varargin{4};
           obj.Imbalanced=varargin{5};
           obj.QTable_node=varargin{6};
           end
           % nobound: The number of boundaries
           % Disthandle: The handle of distance function of each boundary
           % in the order of NoBound
           % Format of BCTable_node:
           % 1st column: The type of the dof, 1- xdisp, 2-ydisp, 3-pressure
           % 2nd column: The global index of the node
           % 3rd and 4th: The coordinates of the node
           % 5th column: The prescribed displacement or pressure value.
           % BCTable_line: The boundary condition type matrix, each row
           % represents one boundary, each column represents one dof
                % [xdisp,ydisp, pressure and loads], MAKE SURE THE LAST COLUMN IS FOR LOADS CONDITION
                %-- For the first three rows
                % '0' means no need to explicitly apply boundary condition;
                % '1' means predefineded nodal value (displacement or pressure)
                % Specially, for the third row
                % '0' also means impermeable condition
                % '5' means volumetric normal flow boundary condition, qn,
                % outward positive
                % '6' means volumetric x-y flow boundary condition [qx,qy]
                %-- For the fourth row
                % '0' means no traction boundary condition;
                % '2' means constant distributing stress loads in cartesian
                % coordinate system, like in-situ stress to the outer boundary
                % '3' means traction boundary condition, [Tnormal,Tshear].
                % Tnomal is related to the stress and the outward unit vector
                % '4' means traction boundary condition, [Tx,Ty].
           % BCTable_node: boudnary condition type matrix of each node with
           % BCTable_disp: structured cell for each boundary, each cell
           % contains two rows, first row is node index, second row the
           % prescribed dof;
           % BCTable_stress: 6 columns, first two rows represent the
           % nodes along the element side in the right order, the rest four
           % rows is the specified stress tensor in vector form.
           % BCTable_traction: 4 columns. first two rows represent the
           % nodes along the element side in the right order, the third row
           % should be the x component and the fourth row should be the
           % y component.
           % BCTable_flow: 3 columns. first two rows represent the
           % nodes along the element side in the right order, the third row
           % should be the normal flow.
           % special care: [idof,id,x-coord,y-coord,disp]

       end
      
       function preprocess(obj)
           if obj.Mesher.Nface==3
           obj.Mesher.validate;
           obj.Mesher.reorder;
           end       
           obj.Mesher.nodetype;
           obj.BoundEdge=boundedges2(obj.Mesher.p,obj.Mesher.t);
           obj.constrBCTable;
           obj.crtpsdbound;
%            save('preprocessor.mat','obj');
       end    
       constrBCTable(obj);
       crtpsdbound(obj);
   end
end
