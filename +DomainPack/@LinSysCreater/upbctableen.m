function upbctableen(obj)
% should be called for every time step after updating the crack
% geometry (now in domain.updatedofarray_enriched,0301)
% Format of BCTableEn:
% 1st column: The index of dof
% 2nd column: The prescribed displacement or pressure value
% BCTableEn=[35227,0;35230,0];                            % start over, do not appendix on the old BCTableEn.
BCTableEn=[];
% this is specific to mesh-elemfile_01132020_carrier.txt, to ensure the
% ux_enr on the left edge also equals to zero. 01212020
for iEnrich=1:length(obj.EnrichItems)
    enrichitem=obj.EnrichItems{iEnrich};
    id=enrichitem.Id;
    stdnodes=enrichitem.Mygeo.Stdnodes;
%     bctable=zeros(length(stdnodes)*3,2);
    bctable=zeros(length(stdnodes)*2,2);
    for istd=1:length(stdnodes)
%         bctable(3*istd-2:3*istd,1)=obj.NodeDict(stdnodes(istd)).EnrDofArray{id};
        bctable(2*istd-1:2*istd,1)=obj.NodeDict(stdnodes(istd)).EnrDofArray{id}(1:2);
    end
    BCTableEn=[BCTableEn;bctable];
    
end
obj.BCTableEn=BCTableEn;
end
