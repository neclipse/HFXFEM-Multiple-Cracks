function crtLHS_UP(obj,newmark,iinc)
% Unify the LHS assembly using the single rearraged elemental StifMatrix;
% The global dofs are arragned as
% [ux1,uy1,p1,ux2,uy2,p2....uxend,uyend,pend,uxenrN,uyenrN,penrN,...],
% which is consistent with elemental elem.Locarray.
% However,the elem.StifMatrix corresponds to dofs
% [ux1,uy1,ux2,..ux4,uy4,p1,...,p4,uxenr1,uyenr1,...uxenr4,uyenr4,penr1,...,penr4];
% Therefore, it is not the original locarry should be used in the assembly
% of LHS. Instead LocarrayU(1:2:end-1)=Locarray(1:3:end);
% LocarrayU(2:2:end)=Locarray(2:3:end); LocarrayP=Locarray(3:3:end);
% Locarray=[LocarrayU,LocarrayP];
Ntriplets=numel(obj.ElemDict)*24*24;       % The maximum number of components of I, J and Value, one node has 6 dofs in average.
I=zeros(Ntriplets,1);
J=zeros(Ntriplets,1);
Value=zeros(Ntriplets,1);
Valuenew=zeros(Ntriplets,1);
ntriplets=0;        % counter
gbinp=obj.ElemDict(1).GaussPntDictM(1).GBINP;
%% Temporarily add this module to differentiate all blending elements from standard elements
blendingelems=zeros(length(obj.ElemDict),1);
icount=1;
for ienr=1:length(obj.EnrichItems)
    elems=obj.EnrichItems{ienr}.Mygeo.Blendingelems;
    icount2=icount+length(elems)-1;
    blendingelems(icount:icount2)=elems;
    icount=icount2+1;
end
blendingelems=blendingelems(1:icount2);
blendingelems=unique(blendingelems);
for ielem=1:length(obj.ElemDict)
    % fprintf('Assembling the %d element\n',ielem);   
    % enriched element
    if any(blendingelems==ielem)
        % blending element
        blending=true;
    else
        % standard element
        blending=false;
    end
    ienrich= find(obj.ElemDict(ielem).Enrich);
    if ~isempty(ienrich)
        % there is at most only one ienrich at current stage 11/22/18
        obj.ElemDict(ielem).crtstif_enriched_elemental(newmark,ienrich,gbinp,blending);
        Std=obj.ElemDict(ielem).Locarray;
        Enr=obj.ElemDict(ielem).JacobianMat.LocarrayEnr;
        stdU=zeros(size(obj.ElemDict(ielem).JacobianMat.LocarrayU));
        enrU=zeros(size(obj.ElemDict(ielem).JacobianMat.LocarrayUEnr));
        stdU(1:2:end-1)=Std(1:3:end);
        stdU(2:2:end)=Std(2:3:end);
        stdP=Std(3:3:end);
        enrU(1:2:end-1)=Enr(1:3:end);
        enrU(2:2:end)=Enr(2:3:end);
        enrP=Enr(3:3:end);
        Locarray=[stdU,stdP,enrU,enrP];
        % A row or a column of WJ follows the order of [stdU, stdP, enrU, enrP]
        Jacobian=obj.ElemDict(ielem).JacobianMat.WJ;
        Jacobiannew=obj.ElemDict(ielem).JacobianMat.WJnew;
    else
        obj.ElemDict(ielem).crtstif(newmark,iinc,gbinp,blending);
        Allgloballocs=obj.ElemDict(ielem).Locarray;
        LocarrayU=zeros(size(obj.ElemDict(ielem).LocarrayU));
        LocarrayU(1:2:end-1)=Allgloballocs(1:3:end);
        LocarrayU(2:2:end)=Allgloballocs(2:3:end);
        LocarrayP=Allgloballocs(3:3:end);
        Locarray=[LocarrayU,LocarrayP];
        % A row or a column of StifMatrix follows the order of [stdU, stdP]
        Jacobian=obj.ElemDict(ielem).StifMatrix;
        Jacobiannew=Jacobian;
    end
%     % Assembly start
    row=Locarray;
    col=row;
%     for irow=1:length(row)
%         for icol=1:length(col)
%             ntriplets = ntriplets + 1 ;
%             I(ntriplets) = row(irow); % row index in global stiffness matrix given by elemental location array
%             J(ntriplets) = col(icol);
%             Value(ntriplets) =Jacobian(irow,icol); % value extracted from the elemental Jacobain matrix
%             Valuenew(ntriplets) =Jacobiannew(irow,icol); % value extracted from the elemental Jacobain matrix
%         end
%     end
    % It is claimed that MATLAB is more efficient in looping by columns
    for icol=1:length(col)
        for irow=1:length(row)
            ntriplets = ntriplets + 1 ;
            I(ntriplets) = row(irow); % row index in global stiffness matrix given by elemental location array
            J(ntriplets) = col(icol);
            Value(ntriplets) =Jacobian(irow,icol); % value extracted from the elemental Jacobain matrix
            Valuenew(ntriplets) =Jacobiannew(irow,icol); % value extracted from the elemental Jacobain matrix
        end
    end
end
% To remove extra zeros
I=I(1:ntriplets);
J=J(1:ntriplets);
Value=Value(1:ntriplets);
Valuenew=Valuenew(1:ntriplets);
JMatrix=sparse(I,J,Value,obj.Drow,obj.Dcol);
JMatrixnew=sparse(I,J,Valuenew,obj.Drow,obj.Dcol);
% To get the guaranteed symmetric matrix (applied on 02102019). 
%(Turn off on 02142019 as of the systematic asymmetry for partially
%enriched elements)
% (03032019, turned on again because the symmetry for the partially
% eneriched elements are fixed)
A=tril(JMatrix);
B=A+tril(A,-1)';
obj.LHS=B;
Anew=tril(JMatrixnew);
Bnew=A+tril(Anew,-1)';
obj.LHSnew=Bnew;

% obj.LHS=JMatrix;
% obj.LHSnew=JMatrixnew;
end

