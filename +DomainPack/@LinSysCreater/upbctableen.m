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
for iEnrich=1:length(obj.EnrichItems)
    enrichitem=obj.EnrichItems{iEnrich};
    id=enrichitem.Id;
    stdnodes=enrichitem.Mygeo.Stdnodes;
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
