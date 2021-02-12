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
       NewNodes         % Newly found nodes to be enriched at this increment
       NewElems         % These info will be used in obj.postprocess for skipping.
       TransElems   % Part of the NewElems that are actually transited from previously smeared element.
       InitialMode      % Initial crack mode: 1-perforated; 2-smeared 3-compressive; 4-tensile.
       Smeared=false;   % Boolean flag to tell if the initial crack is smeared, dependent on InitialMode
       Cohesive         % flag for the type of TSL, 'linear'-linear softening; 'bilinear'-bilinear softening.
       Alpha=pi/2;      % angle between the crack plane and the initial loading for inplace mode
       Isactive = true; % Indicate if the crack is able to propagate from the ti.
   end
   properties(NonCopyable)
       Mytips
   end
   properties(Dependent)
       Stdnodes
   end
   methods
       
       function obj = EnrCrackBody(type,elemdict,nodedict,mygeo,initialmode,cohesive,varargin)
           % To allow constructor work with empty input, 12/07/20.
           if nargin==0
               super_args={};
           else
               super_args={type,elemdict,nodedict};
           end
           obj = obj@EnrichPack.EnrichItem(super_args{:});
           obj.Mygeo=mygeo;
           obj.Id=mygeo.Id; % Mygeo is the universal id of a crack
           obj.Mesh=mygeo.Mesh;
           obj.InitialMode=initialmode;   
           obj.Cohesive=cohesive;       
           obj.generatetips;
           obj.setinteractedelem;
           obj.setenrichednode;
           if ~isempty(varargin)
               obj.Alpha=varargin{1}; % the default value of alpha is pi/2
           end
           % think about how to do with smeared crack...01/31/2021
           % Perhaps, I can set up a derived EnrCrackBody_smeared with almost every
           % attribute of the current EnrCrackBody. With smeared==true, the
           % methods of EnrCrackBody_smeared will be called instead of the
           % not smeared. In this way, it may be easier to transfer
           % attributes between the smeared to open crack. (a proposal) 
           if initialmode==2
               % For now, the newelems and newnodes will be empty.
               obj.Smeared=true;
           else
               % Initially the new elems and nodes are all nodes interacted
               % if the obj.smeared is not true
               obj.NewElems=obj.Interactedelem';
               obj.NewNodes=obj.Enrichednode';
           end
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
       function stdnodes=get.Stdnodes(obj)
           stdnodes=zeros(length(obj.Mytips),2);
           for itip=1:length(obj.Mytips)
               tipelem=obj.Mytips(itip).INTELEM;
               if ~tipelem.Smeared(tipelem.Enrich==obj.Id)
                   stdnodes(itip,:)=obj.Mygeo.Stdnodes(2*itip-1:2*itip);
               end
           end
       end
       
       function checkactive(obj)
           % 12/07/20 The activeness primarily depends on the geometry
           % checking when the package was designed for single crack.
           % It should, however, be updated when there are multiple cracks.
           if isempty(obj.Mygeo.Rtips)
               obj.Isactive = false;
           end
       end
       
       %% function prototyping
       [fe,locarray_enr]=cal_qextenr(obj,q,varargin);
       [unstablegrow,cutflag]=check_grow(obj,varargin);
       initial_enrich_1(obj,varargin);
       initial_enrich_2(obj,varargin);
%        initial_enrich_3(obj,varargin);
       update_enrich_1(obj,varargin);
       update_enrich_2(obj,varargin);
%        update_enrich_3(obj,varargin)
       showme( obj,typex,varargin );
   end
end