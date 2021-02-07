function switching( obj,mode )
%SWITCHING MEthod of Linsyscreater
%   Developed For standard dofs
%   Should not change much for enriched dofs as they are stored together 
%   within the linsyscreater.
%   It turns out only the Cohesive force vector stored in the JacobianMat 
%   of the enriched elements should also be switched between old and new.
if mode==1
    obj.UO=obj.U;                   
    obj.PO=obj.P;
    obj.UOt1=obj.Ut1;
    obj.UOt2=obj.Ut2;
    obj.POt1=obj.Pt1;
    % For the enriched 
    if ~isempty(obj.EnrichItems)
        smeared=[obj.EnrichItems.Smeared];
        realenrich=find(~smeared);
        for i=1:length(realenrich)
            iEnrich=realenrich(i);
            enrichitem=obj.EnrichItems(iEnrich);
            id=enrichitem.Id;
            for ielem=1:length(enrichitem.INTELEM)
                elem=enrichitem.INTELEM(ielem);
                % when the enrichitem.Smeared==false, there is still a
                % chance that some of its involved elements has opened.
                ind=elem.Enrich==id; % should have been fixed for issue #1.
                if ~elem.Smeared(ind) % the crack is open at this element
                    Jacob=elem.JacobianMatDict(ind);
                    % There may be bug when there is change in elem.Enrich.
                    % 02/05/21
                    Jacob.F_coh_old=Jacob.F_coh_i;
                    for ig=1:length(elem.LineGaussDict{ind})
                        elem.LineGaussDict{ind}(ig).TractionO=elem.LineGaussDict{ind}(ig).Traction;
                    end
                end
            end
        end
    end
%% Mode 2 is used for the (dynamic) formulation where the unknowns are incremental dofs, not used by 11/28/18
% elseif mode==2
%     obj.Unknowns=zeros(obj.Drow,1); % NEEDED for upconf
%     obj.U=obj.UO;
%     obj.P=obj.PO;
%     % For the enriched 
%     for iEnrich=1:length(obj.EnrichItems)
%         enrichitem=obj.EnrichItems{iEnrich};
%         id=enrichitem.Id;
%         for ielem=1:length(enrichitem.INTELEM)
%             Jacob=enrichitem.INTELEM(ielem).JacobianMatDict(id);
%             Jacob.F_coh_i=Jacob.F_coh_old;
%         end
%     end
%% Mode 3 is used for the (dynamic) formulation where the unknowns are total dofs, 11/28/18
elseif mode==3
    obj.Unknowns=zeros(obj.Drow,1); % NEEDED for upconf
    obj.U=zeros(obj.Nudof,1);
    obj.P=zeros(obj.Npdof,1);
    if ~isempty(obj.EnrichItems)
        % For the enriched 
        smeared=[obj.EnrichItems.Smeared];
        realenrich=find(~smeared);
        for i=1:length(realenrich)
            iEnrich=realenrich(i);
            enrichitem=obj.EnrichItems(iEnrich);
            id=enrichitem.Id;
            for ielem=1:length(enrichitem.INTELEM)
                elem=enrichitem.INTELEM(ielem);
                % when the enrichitem.Smeared==false, there is still a
                % chance that some of its involved elements has opened.
                ind=elem.Enrich==id;
                if ~elem.Smeared(ind) % the crack is open at this element
                    Jacob=enrichitem.INTELEM(ielem).JacobianMatDict(ind);
                    Jacob.F_coh_i=Jacob.F_coh_old; % initial F_coh_old should be assigned in the JacobianMat.
                    for ig=1:length(elem.LineGaussDict{ind})
                        elem.LineGaussDict{ind}(ig).Traction=elem.LineGaussDict{ind}(ig).TractionO;
                    end
                end
            end
        end
    end
end

end

