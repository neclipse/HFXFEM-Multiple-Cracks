function obj=storage_old( obj, iinc, inc, varargin )
%STORAGE Method of Domain class, to store important results for
%postprocessing.
%   It runs after every converged load increment in the 'running' method. The method will store
%   results into an instance of 'Postprocessor'.
% Prepare an object of postprocessor class

%% NOT UPDATED FOR ENRICHED DOMAIN YET. 03132019 (Finished on 0409)
global inistress inipore;
if nargin==3
    ind=iinc;                                               % index of the post in the postprocessor array is the same as iinc
elseif nargin==4
    ind=varargin{1};                                        % index of the post in the postprocessor array is different from iinc
end
noelem=length(obj.ElemDict);
% gausscol=1;                                               % number of gaussian points to represent one element
post=ToolPack.Postprocessor();
%% Store the nodal results
post.IInc=iinc;
post.Inc=inc;
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
% NOT APPLICABLE TO ENRICHED ELEMENTS 04092019
SN=zeros(2,length(inistress),size(post.XN,1));
ESnN=SN;SPN=SN;
icount=ones(1,length(inistress));
for ielem=1:noelem
    ienrich= find(obj.ElemDict(ielem).Enrich,1); % At most one enrichement in an element at current stage.
    % standard element
    if isempty(ienrich)
        elem=obj.ElemDict(ielem);
        nodes=elem.NodList;
        disp=elem.Un2i;
        p=elem.Pn2i;                                      % Excess pore pressure
        % -- Loop over gaussian points dictionary
        for igauss=1:length(elem.GaussPntDictM);       % number of gaussian points within one element
            elem.GaussPntDictM(igauss)=elem.GaussPntDictM(igauss).matsu(disp,p); % ACCURATE
            %            gauss=elem.GaussPntDictM(igauss);
            %            post.XG(indgauss)=gauss.X;
            %            post.YG(indgauss)=gauss.Y;
            %            gauss=gauss.matsu(disp,p);                 % calculate the strain and total stress values, and total pore pressure at the current gaussian point
            % -- Add the stress values to the corresponding nodes,
            %  note the gaussian points and the nodes for local element follow the same numbering order.
            ESnN(:,:,nodes(igauss))=ESnN(:,:,nodes(igauss))+[elem.GaussPntDictM(igauss).Strainc';icount]; % Use the icount to count the times that to this node the stress is accumulated
            SN(:,:,nodes(igauss))=SN(:,:,nodes(igauss))+[elem.GaussPntDictM(igauss).Stress';icount];
            SPN(:,:,nodes(igauss))=SPN(:,:,nodes(igauss))+[elem.GaussPntDictM(igauss).Stressp';icount];
        end
    else % Enriched elements
        % NOT finished 04092019
        % Need to use enriched interpolation method to determine stress
        % B*disp+ Benr* dispenr
    end
end
% Average the results and reshape them into 2d matrix
ESnN(1,:,:)=ESnN(1,:,:)./ESnN(2,:,:); ESnN(2,:,:)=[];   % averaging the value at the node, NOT ACCURATE
ESnN=permute(ESnN,[1,3,2]);                             % permute the 2nd and 3rd dimension
ESnN=reshape(ESnN,[],length(inistress));                % reshape it to a 2D matrix
SN(1,:,:)=SN(1,:,:)./SN(2,:,:); SN(2,:,:)=[];
SN=permute(SN,[1,3,2]);                                 % permute the 2nd and 3rd dimension
SN=reshape(SN,[],length(inistress));
SPN(1,:,:)=SPN(1,:,:)./SPN(2,:,:); SPN(2,:,:)=[];
SPN=permute(SPN,[1,3,2]);                               % permute the 2nd and 3rd dimension
SPN=reshape(SPN,[],length(inistress));
post.ESnN=ESnN;
post.SN=SN;
post.SPN=SPN;
 
obj.Postprocess(ind)=post;
% Store the practical information for the cracks (0322 not finished, try to use meta.class)
% if ~isempty(obj.EnrichItems)
%    for ienr=1:length(obj.EnrichItems)
%        post.SPN
%    end
% end
end

