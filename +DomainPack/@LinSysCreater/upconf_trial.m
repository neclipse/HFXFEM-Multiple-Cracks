function  upconf_trial( obj,newmark)
% Keep only those terms that explicitly attached to total u-p unknowns 
% at step n+1 at the left side, and move the rest to the R.H.S.
% The operation is equivalent to setting linsys.U and linsys.P to zero and 
% updating their derivatives using the updated newmark coefficients.
% When the inforcer is evaluated and the residual load vector will finish the 
% forementioned operation automatically because here we set Due as zero and
% the f_coh is same as the last increment (see ifstd_enriched).
% Therefore, f_coh evaluated at the end of step n should be directly move to
% the R.H.S because it is not explicitly expressed in terms of total u-p in
% the formulation, K_coh in the Jacobian matrix is only valid for iterative
% duenr (or Duenr for modified newton-raphson scheme) because K_coh is by 
% nature nonlinear. (04152019)
% It is different to treat f_int, q_int and q_intenr because they are
% eventually expressed in total U and P. Noted on 11/30/2018
    nm=newmark;                     % newmark parameters in struct array
%   calculate the trial time derivatives of the updated U and P using
%   newmark method for n+1 time step
%     obj.U=zeros(obj.Nudof,1);       % already done in switching(3)
%     obj.P=zeros(obj.Npdof,1);       % already done in switching(3)
    obj.Ut1=nm.a1*(-obj.UO)-nm.a3*obj.UOt1-nm.a5*obj.UOt2;
    obj.Ut2=nm.a0*(-obj.UO)-nm.a2*obj.UOt1-nm.a4*obj.UOt2;
    obj.Pt1=nm.a1p*(-obj.PO)-nm.a3p*obj.POt1;
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
    % trial update the enriched dofs as well, AT ELEMENTAL LEVEL, 10/12/20
	for iEnrich=1:length(obj.EnrichItems)
		enrichitem=obj.EnrichItems(iEnrich);
% 		id=enrichitem.Id;
		for ielem=1:length(enrichitem.INTELEM)
			Jacob=enrichitem.INTELEM(ielem).JacobianMat;
% 			DofEnrVec=Jacob.LocarrayEnr;
			UEnrVec=Jacob.LocarrayUEnr;
			PEnrVec=Jacob.LocarrayPEnr;
            % For the trial step, the dUn2iEnr should vanish. Because the
            % cohesive force should only corresponds to the crack
            % displacement at the end of last increment 11/29/18
% 			Jacob.dUn2iEnr=zeros(length(UEnrVec),1);    % per interation
            Jacob.DUn2Enr=zeros(length(UEnrVec),1);     % per increment ADDED ON 04032019
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

