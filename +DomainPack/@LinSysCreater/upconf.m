function  upconf( obj,newmark)
%UNPCONF method of linear system creater
%   update configuration: incremental displacements and total accumulated
%   displacements, pressure and all required time derivatives
    ind=zeros(size(obj.Unknowns));
    ind(3:3:end)=1;
    pind=logical(ind);              % logical pore pressure index;
    dispind=~pind;
    du=obj.Unknowns(dispind);       % elemental iterative displacement vector, du_sup(n+1)_sub(i+1)
    dp=obj.Unknowns(pind);          % elemental iterative pore pressure vector, dp_sup(n+1)_sub(i+1)
    % ERROR ON THE UPDATE OF OBJ.U AND OBJ.P, noticed on 6/4, not fixed
    obj.U=obj.U+du;                 % U_sup(n+1)_sub(i+1), displacement
    DU=obj.U-obj.UO;                % incremental change in displacment.
    obj.P=obj.P+dp;                 % p_sup(n+1)_sub(i+1), current excess pore water pressure
    nm=newmark;                     % newmark parameters in struct array
%   calculate the time derivatives of the updated U and P using newmark method
    obj.Ut1=nm.a1*(obj.U-obj.UO)-nm.a3*obj.UOt1-nm.a5*obj.UOt2;
    obj.Ut2=nm.a0*(obj.U-obj.UO)-nm.a2*obj.UOt1-nm.a4*obj.UOt2;
    obj.Pt1=nm.a1p*(obj.P-obj.PO)-nm.a3p*obj.POt1;
%         obj.Ut1=(obj.U-obj.UO)/nm.dt;
%         obj.Ut2=(obj.Ut1-obj.UOt1)/nm.dt;
%         obj.Pt1=(obj.P-obj.PO)/nm.dt;
    nelem=length(obj.ElemDict);
    for ielem=1:nelem
        DofVec=obj.ElemDict(ielem).Locarray;
        UVec=obj.ElemDict(ielem).LocarrayU;
        PVec=obj.ElemDict(ielem).LocarrayP;
        obj.ElemDict(ielem).dXn2i=obj.Unknowns(DofVec);
        obj.ElemDict(ielem).Un2i=obj.U(UVec);
        obj.ElemDict(ielem).Pn2i=obj.P(PVec);   
        obj.ElemDict(ielem).Un2t1=obj.Ut1(UVec);
        obj.ElemDict(ielem).Un2t2=obj.Ut2(UVec);
        obj.ElemDict(ielem).Pn2t1=obj.Pt1(PVec);
    end
    % update the enriched dofs as well, AT ELEMENTAL LEVEL 10/12/20
	for iEnrich=1:length(obj.EnrichItems)
		enrichitem=obj.EnrichItems{iEnrich};
% 		id=enrichitem.Id;
		for ielem=1:length(enrichitem.INTELEM)
			Jacob=enrichitem.INTELEM(ielem).JacobianMat;
% 			DofEnrVec=Jacob.LocarrayEnr;
			UEnrVec=Jacob.LocarrayUEnr;
			PEnrVec=Jacob.LocarrayPEnr;
            Jacob.DUn2Enr=DU(UEnrVec);  % incremental change in uenr, used for cohesion force calculation.
			Jacob.Un2iEnr=obj.U(UEnrVec);
			Jacob.Pn2iEnr=obj.P(PEnrVec);   
			Jacob.Un2t1Enr=obj.Ut1(UEnrVec);
			Jacob.Un2t2Enr=obj.Ut2(UEnrVec);
			Jacob.Pn2t1Enr=obj.Pt1(PEnrVec);
		end
	end

% nnode=length(obj.NodeDict);
% for inode=1:nnode
%    Disparray=obj.NodeDict(inode).DofArray(1:2);     % The potential third index is for the global location of pressure dof
%    obj.NodeDict(inode).ItrDispList=obj.Unknowns(Disparray);
%    obj.NodeDict(inode).IncDispList=obj.NodeDict(inode).IncDispList+obj.NodeDict(inode).ItrDispList;
%    obj.NodeDict(inode).TotDispList=obj.NodeDict(inode).TotDispList+obj.NodeDict(inode).ItrDispList;
% end
end

