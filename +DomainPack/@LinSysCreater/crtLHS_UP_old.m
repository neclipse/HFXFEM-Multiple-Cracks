function crtLHS_UP_old(obj,newmark,iinc)
% This is the old approach as a new approach was developed when developing
% the assembly for enriched system. It turns out this method is not optimal
% when there are too many sub-matrices. 11/25/2018

% IMPORTANT ERROR: WITHOUT THIS CALCULATION, THE FOLLOWING I J
% VALUE WOULD BE EMPTY.

obj.ElemDict(1).crtstif(newmark,iinc);
[dr,~]=size(obj.ElemDict(1).StifMatrix);
[dr11,dc11]=size(obj.ElemDict(1).JM11);
Nntriplets11=length(obj.ElemDict)*dr11*dc11;
I11=zeros(Nntriplets11,1);
J11=I11;Value11=I11;
[dr12,dc12]=size(obj.ElemDict(1).JM12);
Nntriplets12=length(obj.ElemDict)*dr12*dc12;
I12=zeros(Nntriplets12,1);
J12=I12;Value12=I12;
[dr21,dc21]=size(obj.ElemDict(1).JM21);
Nntriplets21=length(obj.ElemDict)*dr21*dc21;
I21=zeros(Nntriplets21,1);
J21=I21;Value21=I21;
[dr22,dc22]=size(obj.ElemDict(1).JM22);
Nntriplets22=length(obj.ElemDict)*dr22*dc22;
I22=zeros(Nntriplets22,1);
J22=I22;Value22=I22;
ntriplets11=0;        % total count of components in J11 the matrix
ntriplets12=0;        % total count of components in J12 the matrix
ntriplets21=0;        % total count of components in J21 the matrix
ntriplets22=0;        % total count of components in J22 the matrix
% Map index from elemental local J to global Jmatrix
[mapug,mappg]=mapping(dr11,dr22,dr);
for ielem=1:length(obj.ElemDict)
    % fprintf('Assembling the %d element\n',ielem);
    obj.ElemDict(ielem).crtstif(newmark,iinc);
    row=obj.ElemDict(ielem).Locarray;
    col=row;
    for irow=1:dr11
        for icol=1:dc11
            ntriplets11 = ntriplets11 + 1 ;
            I11(ntriplets11) = row(mapug(irow)); % row index in global stiffness matrix given by location array
            J11(ntriplets11) = col(mapug(icol));
            Value11(ntriplets11) = obj.ElemDict(ielem).JM11(irow,icol) ;
        end
    end
    for irow=1:dr12
        for icol=1:dc12
            ntriplets12 = ntriplets12 + 1 ;
            I12(ntriplets12) = row(mapug(irow)); % row index in global stiffness matrix given by location array
            J12(ntriplets12) = col(mappg(icol));
            Value12(ntriplets12) = obj.ElemDict(ielem).JM12(irow,icol) ;
        end
    end
    for irow=1:dr21
        for icol=1:dc21
            ntriplets21 = ntriplets21 + 1 ;
            I21(ntriplets21) = row(mappg(irow)); % row index in global stiffness matrix given by location array
            J21(ntriplets21) = col(mapug(icol));
            Value21(ntriplets21) = obj.ElemDict(ielem).JM21(irow,icol) ;
        end
    end
    for irow=1:dr22
        for icol=1:dc22
            ntriplets22 = ntriplets22 + 1 ;
            I22(ntriplets22) = row(mappg(irow)); % row index in global stiffness matrix given by location array
            J22(ntriplets22) = col(mappg(icol));
            Value22(ntriplets22) = obj.ElemDict(ielem).JM22(irow,icol) ;
        end
    end
end
JMat11=sparse(I11,J11,Value11,obj.Drow,obj.Dcol);
JMat12=sparse(I12,J12,Value12,obj.Drow,obj.Dcol);
JMat21=sparse(I21,J21,Value21,obj.Drow,obj.Dcol);
JMat22=sparse(I22,J22,Value22,obj.Drow,obj.Dcol);
JMatrix=JMat11+JMat12+JMat21+JMat22;
obj.LHS=JMatrix;
end
