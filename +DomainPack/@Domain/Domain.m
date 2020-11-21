% This class is called Domain. It plays a role as the boss in the
% builder design pattern.Domain owns everything but does not do everything.
% One object of a class Domain is used to define one problem. 
% LinSysCreater is the CEO of linear system building which supervises the
% worker(class element). Element is the basic builder which produces each
% element of the linear system and then assemle them into a global linear
% system.
classdef Domain
   properties
      Preprocess    % Instance of ToolPack.Preprocessor
      NewtonRaphson % Solver--Newton_Raphson Iterator, the director for elasto-plastic problem
      Postprocess   % Instance array of ToolPack.Postprocessor
      EnrichItems   % Dictionary of the EnrichItem objects
      LinSysCrt
      ElemType      % Specify the element type of the problem;
      % 1-plane stress, 2-plane strain, 3-axisymmetric, 4 UP element
      ElemDict      % Element dictionary, the builder of stiffness matrix
      NodeDict      % Node dictionary, providing the location arrary
      MatType       % Specify the Material type (Constitutive law)
      % 1-linear elastic, 2-Drucker-Prager, 3-Drucker-Prager/Cap, 
      % 4-Coupled material
      NoDofs        % Number of all standard Dofs
      NoUDofs       % Number of all standarddisplacement dofs
      NoPDofs       % Number of all standard pore pressure dofs
      NoEnrDofs=0;     
      NoUenrDofs=0;
      NoPenrDofs=0;
      PsdDofs       % list of all standarddofs on pseudo boundaries
      PsdEnrDofs    % list of all enrddofs detected by the EnrichItems.
      GBINP         % Struct arrary of common used global constants
   end
   
   methods

       function obj=Domain(elemtype,material) % Construct an object of Domain class
           if nargin>0
           obj.ElemType=elemtype;
           obj.MatType=material;
%            obj.GBINP=assembleglobalinputs();
           end
       end
       
       % prototype undefined functions 
       obj=updatedofarray(obj);
       obj=updatedofarray_enriched(obj);
       initiate_enrich(obj);
       obj=assignmesh(obj);     % Use mesh info to initiate elemdict and nodedict
       gaussdictm=crtgpd( obj, varargin );
       obj=crtlinsys(obj);      % create linear equation system
       obj=updatelinsys(obj);
       obj=storage( obj, postdict, iinc, inc, varargin );
       obj=running(obj,postdict,savemode,varargin);
       assign_arbitrary_flow(obj);
       snapshot( obj, variable, timesteps, varargin);
       %
   end
end
