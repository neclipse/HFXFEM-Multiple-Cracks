classdef LinSysCreater < handle
   properties
       LHS              % Left hand side of the linear system
       LHSnew           % Left hand side of the linear system for the new incremnt(designed for cohesive crack 04122019)
       LHSO             % last converged left hand side
       INI              % Initial force vector resulted from initial stress
       RHS              % Right hand side of the linear system
       Drow             % The size of the linear system in row dimension (mixed)
       Dcol             % The size of the linear system in column dimension (mixed)
       Nudof            % (mixed)
       Npdof            % (mixed)
       Unknowns         % The (mixed) unknowns vector, [U_sub(n+1),P_sub(n+1)]
       U                % Current total displacements
       Ut1              % 1st order time derivative of current displacements
       Ut2              % 2nd order
       UO               % Last Converged  displacements
       P                % current excess pore pressure
       Pt1
       PO               % last converged excess pore pressure 
       UOt1             % 1st order time derivative of last converged displacement
       UOt2             % 2st order time derivative of last converged displacement
       POt1             % 1st order time derivative of last converged excess pore pressure
       MatType          % The material type, same as that in Domain class
       ElemType         % The element type, same as that in Domain class
       ElemDict         % The element dictionary created in the domain class
       NodeDict         % The node dictionary created in the domain class
       BCTable_line
       BCTableim        % The imbalanced boundary edges for initial_RHS
       BCTables         % The boundary condition table for stress loads
       BCTablet         % The boundary condition table for traction loads
       BCTablen         % The boundary condition table for sparse nodes
       BCTabled         % The boundary condition table for displacement
       BCTableq         % The boundary condition table for flow loads
       BCTableEn        % The boundary condition table for the sparse standard nodes described in the enrichitems
       BCTableqEn       % The boundary condition table listing the crack mouth inflow for each enrichitem
	   EnrichItems
   end
   methods
       function obj=LinSysCreater(ndof,nudof,npdof,elemtype,mattype,elemdict,...
               nodedict,bctable,bctableim,bctables,bctablet,bctablen,bctabled, ...
			   bctableq,enrichitems)
           obj.Drow=ndof;               % All mixed dofs
           obj.Dcol=ndof;
           obj.Nudof=nudof;             % All displacement dofs including the enriched
           obj.Npdof=npdof;             % All pressure dofs including the enriched
           obj.ElemType=elemtype;
           obj.MatType=mattype;
           obj.ElemDict=elemdict;       % Handles of all elements
           obj.NodeDict=nodedict;       % Handles of all nodes
           obj.BCTable_line=bctable;
           obj.BCTableim=bctableim;
           obj.BCTables=bctables; 
           obj.BCTablet=bctablet; 
           obj.BCTablen=bctablen;
           obj.BCTabled=bctabled;
           obj.BCTableq=bctableq;
           obj.UO=zeros(obj.Nudof,1);
           obj.UOt1=zeros(obj.Nudof,1); 
           obj.UOt2=zeros(obj.Nudof,1);
           obj.PO=zeros(obj.Npdof,1);
           obj.POt1=zeros(obj.Npdof,1);
           obj.EnrichItems=enrichitems;	% Enrichement items should be iterated to assemble the enriched dofs to the whole global system
       end
       %%-- Declare other methods which are not defined in the main file
       solve(obj,newmark,varargin);              % solve the linear system
       crtLHS_UP(obj,newmak,iinc);      % creat the left-hand side of the linear system by assembling the element stiffness matrix
       crtRHS_UP(obj,newmark)           
       appbound_v4(obj);                % apply boundary condition to the linear system
       initialRHS( obj );
       initialize( obj );               % not used in the current version 03012019
       switching(obj,mode);
       upconf(obj,newmark);             % update configuration within one increment
       upconf_trial(obj,newmark);       % update configuration for next increment
       % Disenrich the standard nodes detected by the Enrichitems (03012019). or prescribe the nonzero enriched values (later)   
       upbctableen(obj);
       calelemstress(obj);
   end

end