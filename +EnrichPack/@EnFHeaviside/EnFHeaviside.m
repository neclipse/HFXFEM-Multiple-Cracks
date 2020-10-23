classdef EnFHeaviside < EnrichPack.EnrichFun
   properties
       Type                                 % type of heaviside function, 
       % 1. Standard heaviside, 2.Sign-step function
       % Different type of heaviside function will lead to a different
       % formulation of displacement jump. 
       % The first type renders N*u_tilde while
       % the second type renders 2*N*u_tilde. However, the magnitude of
       % t_tilde will automatically adjust itself which lead to a
       % equivalent displacement jump. To avoid the confusing definition,
       % it is by default adopting the first type of heaviside function.
       % 05072019, for the second type, divide the sign-step by 2, then the
       % formulation for displacement jump will be the same: Nstd*u_enr
       Smoothed                             % flag true means, the function is smoothed
       Epsilon                              % tiny number around zero with respect to element size
       Lsv                                  % lsv is Mygeo.Phi: (nodes,phi)
   end
   
   methods
       function obj=EnFHeaviside(smoothed,minelength,lsv,varargin)
           if isempty(varargin)
               obj.Type=1;
           else 
               obj.Type=varargin{1};
           end
           obj.Smoothed=smoothed;         
           %% IMPORTANT BUG: OBJ.EPSILON IS A CRITICAL PARAMETER TO THE XFEM
           % FRAMEWORK, IT SHOULD BE DETERMINED BASED ON THE MESH AND THE
           % SUBDOMAIN METHOD SO THAT ONLY THE CLOSEST GAUSSIAN POINTS TO
           % THE CRACK ARE WITHIN [-OBJ.EPSILON,+OBJ.EPSILON] IN THE
           % DETERMINATION OF DERIVATIVES. 01272020.
           %% Confirmed that the obj.Epsilon should be very small so that 
           % NO gaussian points are included in the range. 013120
           obj.Epsilon=0.0001*minelength;    % minelength is the minimun element length along the crack   
           obj.Lsv=lsv;
       end
       
       function enf=calculate(obj,phi)
           % purely calculate the enrichment function value(s) based on the
           % given phi value(s)
           switch obj.Type
               case 1 %[-1,1]
                   if ~obj.Smoothed
                       enf=heaviside(phi);
                   else
                       enf=zeros(size(phi));
                       L1=phi<-obj.Epsilon;
                       L2=phi>obj.Epsilon;
                       L3=~L1 & ~L2;
                       enf(L2)=1;
                       enf(L3)=1/2+phi(L3)/2/obj.Epsilon+1/2/pi*sin(pi*phi(L3)/obj.Epsilon);
                   end
               case 2 % [-1/2, 1/2]
                   if ~obj.Smoothed
                       enf=sign(phi)/2; 
                       % it is important to be divided by 2 otherwise, the
                       % crackdisplacement expression becomes 2*Nu*uenr,
                       % which may be used in multiple places like
                       % crtstif_enriched and updateopening.
                   else
                       enf=ones(size(phi));
                       L1=phi<-obj.Epsilon;
                       L2=phi>obj.Epsilon;
                       L3=~L1 & ~L2;
                       enf(L1)=-1/2;
                       enf(L3)=(phi(L3)/obj.Epsilon+1/2/pi*sin(pi*phi(L3)/obj.Epsilon))/2;                                                                  
                   end
               otherwise
                   error('The specified Heaviside function type is not recognized')
           end
       end
       function enf_der=calculate_der(obj,phi)
           % purely calculate the enrichment function value(s) based on the
           % given phi value(s)
           switch obj.Type
               case 1
                   if ~obj.Smoothed
                       enf_der=dirac(phi);
                   else
                       enf_der=zeros(size(phi));
                       L1=phi<-obj.Epsilon;
                       L2=phi>obj.Epsilon;
                       L3=~L1 & ~L2;
                       enf_der(L3)=(1+cos(pi*phi(L3)/obj.Epsilon))/2/obj.Epsilon;
                   end
               case 2
                   if ~obj.Smoothed
                       enf_der=dirac(phi);
                   else
                       enf_der=zeros(size(phi));
                       L1=phi<-obj.Epsilon;
                       L2=phi>obj.Epsilon;
                       L3=~L1 & ~L2;
                       enf_der(L3)=3*(1-(phi(L3)/obj.Epsilon).^2)/4/obj.Epsilon;                                                                  
                   end
               otherwise
                   error('The specified Heaviside function type is not recognized')
           end
       end
       
       function enrichelem(obj,elem,id,varargin)
           %% Enrich the elem for the id crack object. The core is obj.enrichgauss
           % can provide the nodes_phi directly
           if ~isempty(varargin)
               nodes_phi=varargin{1};
           else
               % obtain the nodes_phi from obj.Lsv.
               % enrich the gauss points inside
               nodes=elem.NodList;
               %            stdnodes=elem.NodstdList;
               %            stdnodes=stdnodes{id};                  % standard nodes list for the current enrichitem id;
               nodespool=obj.Lsv(:,1);                  % all nodes in the interacted elems
               phipool=obj.Lsv(:,2);
               [Lia,Lcob]=ismember(nodes,nodespool);    % find the relative index of nodes in the nodespool
               if ~all(Lia)
                   error('The requested nodes are not stored in the nodespool');
               end
               nodes_phi=phipool(Lcob);
           end
           % loop over all Gauss points for the current enrcrack(id)
           k=elem.Enrich==id;
           GaussPnt_domain=obj.enrichgauss(nodes_phi,elem.EnrichGauss,k);
           GaussPnt_line=obj.enrichgauss(nodes_phi,elem.LineGaussDict{k});
           % Because gausspnt is an data object, hard copy is required
           elem.EnrichGauss=GaussPnt_domain;
           elem.LineGaussDict{k}=GaussPnt_line;
       end
       
       function addnodedofs(~,node,varargin)
           % enrich the node object by adding corresponding dof and
           % dofarray, not necessary after 10/02/20
           if ~isempty(varargin)
              upindicator=varargin{1};
           else
               upindicator=1;           % u indicator
           end
           if upindicator==1
               node.NoUenrDofs=2;
           elseif upindicator==2
               node.NoPenrDofs=1;
           end
           node.NoEnrDofs=node.NoUenrDofs+node.NoPenrDofs;
       end
       
       function [GaussPnt,Nuenrplus,Nuenrminus]=enrichgauss(obj,nodes_phi,GaussPnt,varargin)
           if isempty(varargin)
               mode=1; % line gaussian points
           else
               mode=2;  % domain gaussian points
               k=varargin{1}; % k is the logical index for the current EnrichItem
           end
          for igauss=1:length(GaussPnt)
               N=GaussPnt(igauss).Np;
               Nu=GaussPnt(igauss).Nu;
               Bmat=GaussPnt(igauss).Bmat;
               % ADDED ON 11/05/18 TO ENSURE LINE GAUSSIAN POINTS HAVE ZERO
               % GAUSS_PHI. NO NEED TO CHANGE RIDGE FUNCTION AS THAT IS
               % ALREADY CONTINUED.
%                if igauss<(length(GaussPnt)-nline+1)
                   gauss_phi=N*nodes_phi;
%                else
%                    gauss_phi=0;
%                end
               % shifted Phi for every enriched node
               heaviside_enrnodes=obj.calculate(nodes_phi);
               heaviside_gauss=obj.calculate(gauss_phi);
               Phishift=heaviside_gauss-heaviside_enrnodes;
               phi_der=obj.calculate_der(gauss_phi); %03282019
               %% On 02262019, new strategy does not delete the standard
               % nodes in this stage, but enforce the enriched dofs as
               % zero in the boundary conditions.
               % delete the standard nodes in the enrichment, IMPORTANT STEP TO ONLY ENRICH THE
               % ENRICHED NODES
%                Phishift(stdnodes)=0;
               if isempty(GaussPnt(igauss).Enf{k})
                   enf=struct('Phi',nan,'Phishift',nan,'Ridge',Ridge,'JunctionU',nan,'JunctionP',nan);
                   GaussPnt(igauss).Enf{k}=enf;
               else
                   GaussPnt(igauss).Enf{k}.Phi=gauss_phi;
                   GaussPnt(igauss).Enf{k}.Phishift=Phishift;
               end
               % prepare Nuenr for all; and Nuenrplus, minus for linegauss
               Nuenr=zeros(size(Nu));
               nuenr=N.*Phishift';
               Nuenr(1,1:2:end)=nuenr;
               Nuenr(2,2:2:end)=nuenr;
               k=find(k); % return the index of logical "1"
               startind=1+(k-1)*size(Nuenr,2);
               endind=k*size(Nuenr,2);
               GaussPnt(igauss).Nuenr(:,startind:endind)=Nuenr;
               % if gausspnt is cohesive class, then assign plus and minus
               % NOTE PLUS AND MINUS ARE ONLY CALCULATED FOR THE CURRENT
               % CRACK ALTHOUGH THEY TAKE THE SAME SHAPE OF NUENR. IT IS
               % JUST TO BE COMPATIBLE TO THE CALL OF NUENRPLUS AND UENR
               % LATER ON IN GAUSS.MATCTU AND GAUSS.UPDATECRACKOPENING.
               if mode==1
                   % To calculate the displacement right above and below
                   % the crack
                   Phishiftplus=obj.calculate(1)-heaviside_enrnodes;
                   Phishiftminus=obj.calculate(-1)-heaviside_enrnodes;
                   Nuenrplus=Nuenr;
                   Nuenrminus=Nuenr;
                   Nuenrplus(1,1:2:end)=N.*Phishiftplus';
                   Nuenrplus(2,2:2:end)=N.*Phishiftplus';
                   Nuenrminus(1,1:2:end)=N.*Phishiftminus';
                   Nuenrminus(2,2:2:end)=N.*Phishiftminus';
                   GaussPnt(igauss).Nuenrplus(:,startind:endind)=Nuenrplus;
                   GaussPnt(igauss).Nuenrminus(:,startind:endind)=Nuenrminus;
               end
               % prepare Bmatenr
               Bmatenr=zeros(size(Bmat));
               N_x=Bmat(1,1:2:end);
               N_y=Bmat(2,2:2:end);
               % The dirac delta function is usually zero because the gauss
               % integration points are designed to be away from the
               % crack.03282019.
               %% IMPORTANT BUG: THE PHISHIFT_XDER OR YDER PART SHOULD NOT 
               %MANUALLY DROPPED, IT MAY BE MASKED WHEN GAUSSIAN POINTS ARE
               %TRULY FAR AWAY FROM THE DISCONTINUITY, BUT IT MAY HAVE AN
               %IMPACT WHEN GAUSSIAN POINTS ARE CLOSE TO THE DISCONTINUITY.
               %012420. However, the gaussian points should not be too
               %close to the discontinuity.
               phishift_xder=phi_der*N_x*nodes_phi;
               phishift_yder=phi_der*N_y*nodes_phi;
               N_xenr=N_x.*Phishift'+N*phishift_xder;
               N_yenr=N_y.*Phishift'+N*phishift_yder;
               Bmatenr(1,1:2:end)=N_xenr;
               Bmatenr(2,2:2:end)=N_yenr;
               Bmatenr(3,1:2:end)=N_yenr;
               Bmatenr(3,2:2:end)=N_xenr;
               startind=1+(k-1)*size(Bmatenr,2);
               endind=k*size(Bmatenr,2);
               GaussPnt(igauss).Bmatenr(:,startind:endind)=Bmatenr;
           end 
       end
   end
end
