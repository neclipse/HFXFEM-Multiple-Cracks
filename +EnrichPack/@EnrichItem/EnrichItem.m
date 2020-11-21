classdef EnrichItem < matlab.mixin.Copyable
   properties
       Type                                                                 % crack body, crack tip, inclusion, void, etc
       Id                                                                   % Id of the same type EnrichItem, automatically updated when new item is constructed
       Myname                                                               % an optional name of the enrichment item, like "left edge crack", or "center void", etc                                              
       Interactedelem
       Enrichednode
%        Allnodes 
       Timestep
   end
   properties(NonCopyable)
       Mesh
       Elemdict
       Nodedict                                                               % vector, a series of discrete time step, indicating the time line of the evolution of enrichment
       Mygeo
       Myenfs
       Myassigner                                                             % all nodes within the enriched elements including the standard nodes
       INTELEM                                          
       ENNODE
   end
   properties(Dependent)
       Allnodes
   end
   
   methods(Abstract)
       setinteractedelem(obj)                                              % Find the interacted elements and assign to the obj.InteractedElem, implemented by mygeo
       setenrichednode(obj)                                                % Select proper nodes upon the knowledge of interacted elements to enrich, implemented by assigner and myenfs
%        initial_enrich(obj)                                                 % initialize the enrichment for this EnrichItem by calling giveinteractedelem, setenrichednode in sequence
%        update_enrich(obj)                                                  % update geometry and associated elements and nodes, selectively enrich new elements and nodes
           % NEED TO MAKE SURE THE ADDED NODES ARE APPENDED AT THE END
           % AFTER CRACK UPDATE
   end

   
   methods
       plotcrack(obj,varargin);
       % constructor
       function obj = EnrichItem(type,elemdict,nodedict) 
%           obj.Id=ID; % assigned by mygeo
          obj.Type  = type;                                                 % type is a string, 'crackbody', 'cracktip','inclusion', 'void',etc
          obj.Elemdict= elemdict;
          obj.Nodedict= nodedict; 
       end
       function set.Mygeo(obj,val)
           if isa(val,'ToolPack.Geometry')
               obj.Mygeo=val;
           else
               error('input must be an instance of ToolPack.Geometry')
           end
       end
%        function set.Mesh(obj,val)
%            if isa(val,'ToolPack.Mesher2d')
%                obj.Mesh=val;
%            else
%                error('input must be an instance of ToolPack.Mesher2d')
%            end
%        end

%      default set, get method. DO not get confused with manual set method
%      such as setinteractedelem
       function set.Myenfs(obj,val)
           for ienf=1:length(val)
               if ~isa(val{ienf},'EnrichPack.EnrichFun')
                   error('input must be an instance of EnrichPack.EnrichFun')
               end
           end
           obj.Myenfs=val;
       end
       function set.Myassigner(obj,val)
           if isa(val,'EnrichPack.EnrichAssigner')
               obj.Myassigner=val;
           else
               error('input must be an instance of EnrichPack.EnrichAssigner')
           end
       end
       function Allnodes=get.Allnodes(obj)
           Allnodes=obj.Enrichednode;
       end
           
       function [totnenrdof, totuenrdof, totpenrdof]=updatedofarray(obj,totnenrdof, totuenrdof, totpenrdof)
           % The initial input of totnenrdof is the number of all standard
           % dofs, i.e., the enriched dofs are indexed after all standard
           % dofs.
  
           % THE ORDER MAY MATTER BECAUSE THE ENRICHED JACOBIAN MATRIX IS 
           % NOT UPDATED BUT THE LOCARRAY MAY CHANGE 11/16/2018
           % NEED TO MAKE SURE THE ADDED NODES ARE APPENDED AT THE END
           % AFTER CRACK UPDATE
           for iN=1:length(obj.NewNodes)
               node=obj.Nodedict(obj.NewNodes(iN));
               % obtain the index where id crack is stored in this
               % element.10/02/20
               k=node.Enrich==obj.Id;
               node.EnrDofArray{k}=totnenrdof+1:totnenrdof+node.NoDofs;
               node.UenrArray{k}=totuenrdof+1:totuenrdof+node.NoUDofs;
               node.PenrArray{k}=totpenrdof+1:totpenrdof+node.NoPDofs;
               % Relieve node.NoEnrDofs and so on by NoDofs because this is
               % not crack-dependent. 10/02/20
               totnenrdof=totnenrdof+node.NoDofs;
               totuenrdof=totuenrdof+node.NoUDofs;
               totpenrdof=totpenrdof+node.NoPDofs;
           end
       end
   end
end