classdef EnrCrackBody < EnrichPack.EnrichItem
   properties
       IntPoints        % Global coordinates of the line gauss points along the whole crack
       Uplus            % Displacement of the upper side u+
       Uminus           % Displacement of the lower side u-
       Aperture         % Apeature at the crack IntPoint
       Pfrack           % Fracture pressures at each the crack IntPoint
       CTraction        % Cohesion traction vector at each crack IntPoint
       CrackVf          % Crack volume fraction on the IntPoint
       CrackVolume      % Integrated crack volume,m^3
       InjectionFlux=0;    % Current injection volumetric rate, m^3/2
       LeakoffFlux=0;      % Current leakoff flux volumetric rate,m^3/2
       LeakoffVolume=0; % Integrated leak volume over time
       Qtable           % Applied fluid injection to this crack
       Length           % Current crack length
       CMOD             % Crack mouth opening
       CMP              % Crack mouth pressure
       % IN such fashion [encrack.Id,q,locnode1,localnode2,x_coord of the
       % injection point,y_coord of the injection point];
       NewNodes
       NewElems
       Perforated       % flag to denote if the existing crack is perforated (completely cohesionless)
       Cohesive         % flag for the type of TSL, 'linear'-linear softening; 'bilinear'-bilinear softening.
       Alpha            % angle between the crack plane and the initial loading for inplace mode
       Isactive = true;
   end
   properties(NonCopyable)
       Mytips
   end
   methods
       
       function obj = EnrCrackBody(id,type,elemdict,nodedict,mygeo,perforated,cohesive,alpha)
           obj = obj@EnrichPack.EnrichItem(id,type,elemdict,nodedict);
           obj.Mygeo=mygeo;
           obj.Mesh=mygeo.Mesh;
           obj.Perforated=perforated;   
           obj.Cohesive=cohesive;       
           obj.Alpha=alpha;
           obj.setinteractedelem;
           obj.setenrichednode;
           obj.generatetips;
           % Initially the new elems and nodes are all nodes interacted
           obj.NewElems=obj.Interactedelem;
           obj.NewNodes=obj.Enrichednode;
           % Later on the new elems and nodes are updated in
           % obj.update_enrich per crack propagation.
       end
       
       function setinteractedelem(obj)                                     % Find the interacted elements and assign to the obj.InteractedElem, implemented by mygeo
            obj.Interactedelem=obj.Mygeo.Intelements;
            obj.INTELEM=obj.Elemdict(obj.Interactedelem);
       end
       
       function setenrichednode(obj)                                       
            obj.Enrichednode=obj.Mygeo.Nodes;
%             obj.Allnodes=[obj.Mygeo.Bodynodes;obj.Mygeo.Rtipnodes];
            obj.ENNODE=obj.Nodedict(obj.Enrichednode);
       end
       
       function set.Mytips(obj,val)
           for itip=1:length(val)
               if ~isa(val(itip),'EnrichPack.EnrCrackTip')
                   error('input must be an instance of EnrichPack.EnrCrackTip')
               end
           end
           obj.Mytips=val;
       end
       
       function initial_enrich(obj)                                         % initialize the enrichment for this EnrichItem
           % assign basic enrichment info to the selected elements and
           % nodes
           elems=obj.Interactedelem;
           nodes=obj.Enrichednode;
           %% set the enrich flag = true and enrich_id=Enrichitem.Id, first
           % nodes then elems
           for iN=1:length(nodes)
               obj.Nodedict(nodes(iN)).setenrich(obj.Id);
           end
           for iE=1:length(elems)
               % set enrich flag and find the standard nodes within the
               % element
               obj.Elemdict(elems(iE)).setenrich(obj.Id);   
               % divide the element into triangular subdomains for 2d
               % integral
               obj.Elemdict(elems(iE)).subdomain(obj.Id);
               % find the gaussian points on the crack for line integral,
               % p=3 to make the gauss quadrature accurate enough for the
               % line integral. (not sure if it is really useful, 03122019)
               % Change p=2 (defautl value) on 06072019.
               obj.Elemdict(elems(iE)).linegauss(obj.Id,obj.Cohesive,obj.Perforated,obj.Alpha);
           end
           %% initial enrich all nodes inside the enriched elements
           % no need to have these after 10/02/20 as udof==2, pdof==1,
           % dofs==3, always.
           % the standard nodes have the enriched dofs but keep zero
%            for iN=1:length(obj.Enrichednode)
%               node=obj.Enrichednode(iN);
%               for ienf=1:length(obj.Myenfs)
%                   myenf=obj.Myenfs{ienf};
%                   myenf.addnodedofs(obj.Nodedict(node),obj.Id);
%               end
%            end
           %% initial enrich elemdict
           for iE=1:length(elems)
               elem=elems(iE);
               % use geometeric info to divide the elem into subdomains
               % for integral purpose.
               %mygeo.subdomain(elem);
               for ienf=1:length(obj.Myenfs)
                   myenf=obj.Myenfs{ienf};
                   myenf.enrichelem(obj.Elemdict(elem),obj.Id);
               end
           end
       end
       
       function checkactive(obj)
           if isempty(obj.Mygeo.Rtips)
               obj.Isactive = false;
           end
       end
       % function prototyping
       
       [fe,locarray_enr]=cal_qextenr(obj,q,varargin);
       [unstablegrow,cutflag]=update_enrich(obj,varargin);
       showme( obj,typex,varargin );
   end
end