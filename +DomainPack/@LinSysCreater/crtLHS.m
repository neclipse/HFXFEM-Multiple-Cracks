function crtLHS(obj,newmark,iinc)
% IMPORTANT ERROR: WITHOUT THIS CALCULATION, THE FOLLOWING I J
% VALUE WOULD BE EMPTY.
obj.ElemDict(1).crtstif(newmark,iinc);
[dr,dc]=size(obj.ElemDict(1).StifMatrix);
Nntriplets=length(obj.ElemDict)*dr*dc;
I=zeros(Nntriplets,1);
J=I;Value=I;
ntriplets=0;        % total count of components in the matrix
for ielem=1:length(obj.ElemDict)
    % fprintf('Assembling the %d element\n',ielem);
    obj.ElemDict(ielem).crtstif(newmark,iinc);
    row=obj.ElemDict(ielem).Locarray;
    col=row; % symmetric square matrix
    for irow=1:dr
        for icol=1:dc
            ntriplets = ntriplets + 1 ;
            I(ntriplets) = row(irow); % row index in global stiffness matrix given by location array
            J(ntriplets) = col(icol);
            Value(ntriplets) = obj.ElemDict(ielem).StifMatrix(irow,icol) ;
        end
    end
end
obj.LHS=sparse(I,J,Value,obj.Drow,obj.Dcol);
end
