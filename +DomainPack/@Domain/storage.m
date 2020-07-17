function obj=storage( obj, postdict, iinc, inc, varargin )
%STORAGE Method of Domain class, to store important results for
%postprocessing.
%   It runs after every converged load increment in the 'running' method. The method will store
%   results into an instance of 'Postprocessor'.
% Prepare an object of postprocessor class

%% NOT UPDATED FOR ENRICHED DOMAIN YET. 03132019 (Finished on 0409)
agl=obj.GBINP;
inipore=agl.inipore;
if nargin==4
    ind=iinc;                                               % index of the post in the postprocessor array is the same as iinc
elseif nargin==5
    ind=varargin{1};                                        % index of the post in the postprocessor array is different from iinc
end
% gausscol=1;                                               % number of gaussian points to represent one element
post=postdict(ind);
%% Store the nodal results
post.IInc=iinc;                                             % Index of the Increment
post.Inc=inc;                                               % The actual time increment value
% Store the reaction force vector on every nodes
a1=obj.NewtonRaphson.Newmark.a1;
% Important to note that the newraphson.IntLoadVect has a coeffient of (-a1)
RLOAD=obj.NewtonRaphson.IntLoadVec(1:obj.NoDofs)/(-a1);
post.RLOADX=RLOAD(1:3:end);
post.RLOADY=RLOAD(2:3:end);
post.XN=obj.Preprocess.Mesher.VX;
post.YN=obj.Preprocess.Mesher.VY;
% Changed on 0409 to adapt to the introduced enrichemnt. Note the unknowns
% on the nodes are directly the solved dofs even for enriched elements
% because enrichment function values on the nodes decrease to zero.
post.UX=obj.LinSysCrt.UO(1:2:obj.NoUDofs);
post.UY=obj.LinSysCrt.UO(2:2:obj.NoUDofs);
post.P=obj.LinSysCrt.PO(1:obj.NoPDofs)+inipore;                            % Total pore pressure = Excess pore pressure+ initial pore pressure
%% Store the crack information 04092019
% Use meta class, actually matlab.mixin.Copyable for automatic clone
enrichitems=cell(size(obj.EnrichItems));
for ien=1:numel(enrichitems)
   enrichitems{ien}=copy(obj.EnrichItems{ien}); 
end
post.EnrichItems=enrichitems;
%% Avearge the integration results to nodal values, NOT ACCURATE, ESPICICALLY ON THE BOUNDARIES
obj.LinSysCrt.calelemstress;
%% Reorder the results
SNtemp=[obj.NodeDict.Stress];
SNPtemp=[obj.NodeDict.Stressp];
SNtemp1=SNtemp(2:5,:);
SNtemp2=repmat(SNtemp(1,:),4,1);
SN=SNtemp1./SNtemp2; 
SNPtemp1=SNPtemp(2:5,:);
SNPtemp2=repmat(SNPtemp(1,:),4,1);
SNP=SNPtemp1./SNPtemp2;
% By old convention, post.SN is length(node)*length(stress)
post.SN=SN';
post.SPN=SNP';
% obj.Postprocess(ind)=post;

end

