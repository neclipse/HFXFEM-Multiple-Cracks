classdef EnFRidge < EnrichPack.EnrichFun
   properties
       Lsv                                  % lsv is Mygeo.Phi: (nodes,phi)
   end
   
   methods
       function obj=EnFRidge(lsv)
           obj.Lsv=lsv;
       end
       
       function enf=calculate(~,nodesphi,N)
           % purely calculate the enrichment function value(s) based on the
           % given nodes phi value(s) and shape functions N
           enf=N*abs(nodesphi)-abs(N*nodesphi);
       end
       function enrichelem(obj,elem,id,varargin)
           if ~isempty(varargin)
               nodes_phi=varargin{1};
           else
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
           GaussPnt_domain=obj.enrichgauss(nodes_phi,elem.EnrichGaussDict{id});
           GaussPnt_line=obj.enrichgauss(nodes_phi,elem.LineGaussDict{id});
           % Because gausspnt is an data object, hard copy is required
           elem.EnrichGaussDict{id}=GaussPnt_domain;
           elem.LineGaussDict{id}=GaussPnt_line;
       end
       function addnodedofs(~,node,id,varargin)
           % enrich the node object by adding corresponding dof and
           % dofarray
           if ~isempty(varargin)
              upindicator=varargin{1};
           else
               upindicator=2;           % p indicator
           end
           if upindicator==1
               node.NoUenrDofs(id)=2;
           elseif upindicator==2
               node.NoPenrDofs(id)=1;
           end
           node.NoEnrDofs(id)=node.NoUenrDofs(id)+node.NoPenrDofs(id);
       end
       function GaussPnt=enrichgauss(obj,nodes_phi,GaussPnt)
           for igauss=1:length(GaussPnt)
               N=GaussPnt(igauss).Np;
               DNp=GaussPnt(igauss).DNp;
               % scalar
               enf=obj.calculate(nodes_phi,N);
               % do not need to use shifted ridge function as ridge_nodes=0
               Ridge=ones(1,GaussPnt(igauss).NNodes)*enf;
               %% On 02262019, new strategy does not delete the standard
               % nodes in this stage, but enforce the enriched dofs as
               % zero in the boundary conditions.
%                Ridge(stdnodes)=0;
               if isempty(GaussPnt(igauss).Enf)
                   enf=struct('Phi',nan,'Phishift',nan,'Ridge',Ridge);
                   GaussPnt(igauss).Enf=enf;
               else
                   GaussPnt(igauss).Enf.Ridge=Ridge;
               end
               % prepare Npenr
               Npenr=N.*Ridge;
               GaussPnt(igauss).Npenr=Npenr;
               % prepare DNpenr
               N_x=DNp(1,:);
               N_y=DNp(2,:);
               % Ridge_x is the partial(Ridge(x))/partial(x),scalar
               Ridge_x=N_x*abs(nodes_phi)-sign(N*nodes_phi)*(N_x*nodes_phi);
               % Ridge_y is the partial(Ridge(x))/partial(x)k,scalar
               Ridge_y=N_y*abs(nodes_phi)-sign(N*nodes_phi)*(N_y*nodes_phi);
               Nj=N;
%                Nj(stdnodes)=0;
               DNpjx=N_x.*Ridge+Nj*Ridge_x;
               DNpjy=N_y.*Ridge+Nj*Ridge_y;
               DNpenr=[DNpjx;DNpjy];
               GaussPnt(igauss).DNpenr=DNpenr;
           end
       end
   end
end