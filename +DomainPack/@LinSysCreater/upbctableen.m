function upbctableen(obj,BCTableqEn)
%% Ensure the ux_enr on the tip edge nodes(stdnodes) equals to zero. 
% EnrDofArray{id}(1:2)are the displacement dofs of the stdnodes.
% should be called for every time step after updating the crack
% geometry (now in domain.updatedofarray_enriched,0301)
% Format of BCTableEn:
% 1st column: The index of dof
% 2nd column: The prescribed displacement or pressure value
% start over, do not appendix on the old BCTableEn.
obj.BCTableqEn=BCTableqEn;
BCTableEn=[];
smeared=[obj.EnrichItems.Smeared];
realenrich=find(~smeared);
for i=1:length(realenrich)
    iEnrich=realenrich(i);
    enrichitem=obj.EnrichItems(iEnrich);
    %     02/05/2021
    % Here is potential bug for the current smeared crack approach: the
    % EnrCrackBody.Smeared has to be consistent with involved elements.
    % Temporary solution: If there is one element in the INTELEM turn open
    % from "smeared", the EnrCrackBody.Smeared will become fale, but the
    % rest of elements.smeared stay true.
    
    % If some of the elements do open later on, the corresponding stdnodes
    % should also be updated. 
    id=enrichitem.Id;
    stdnodes=enrichitem.Stdnodes(:);
    stdnodes=stdnodes(stdnodes~=0); % Now the stdnodes are updated based on the smeared flag 02112021
    %     bctable=zeros(length(stdnodes)*3,2);
    bctable=zeros(length(stdnodes)*2,2);
    for istd=1:length(stdnodes)
        %         bctable(3*istd-2:3*istd,1)=obj.NodeDict(stdnodes(istd)).EnrDofArray{id};
        ind=id==obj.NodeDict(stdnodes(istd)).Enrich;
        bctable(2*istd-1:2*istd,1)=obj.NodeDict(stdnodes(istd)).EnrDofArray{ind}(1:2);
    end
    BCTableEn=[BCTableEn;bctable];
end
obj.BCTableEn=BCTableEn;
end
