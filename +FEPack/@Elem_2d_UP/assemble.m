function [ Gstif ] = assemble( obj,Gstif )
%ASSEMBLE Assemble elemental stiffness matrix into a global stiffness
%matrix with global row and column information calculated by the
%givelocarray function, this is a method called by fullsyscreater. Not used
%by sparsesyscreater.
% row: global row index of dofs 
% col: global col index of dofs
row=obj.Locarray;
col=row;
Gstif(row,col)=Gstif(row,col)+obj.StifMatrix;
end

