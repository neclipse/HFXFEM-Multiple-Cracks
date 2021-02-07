classdef SmearedCrackBody < EnrichPack.EnrCrackBody
   properties
       
   end
   properties(NonCopyable)
       Mytips
   end
   methods
       
       function obj = SmearedCrackBody(type,elemdict,nodedict,mygeo,initialmode,cohesive,varargin)
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
               obj.Smeared=true;
           end
           obj.generatetips;
           obj.setinteractedelem;
           obj.setenrichednode;
           % Initially the new elems and nodes are all nodes interacted
           obj.NewElems=obj.Interactedelem';
           obj.NewNodes=obj.Enrichednode';
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
       initial_enrich_3(obj,varargin);
       update_enrich_1(obj,varargin);
       update_enrich_2(obj,varargin);
       update_enrich_3(obj,varargin)
       showme( obj,typex,varargin );
   end
end