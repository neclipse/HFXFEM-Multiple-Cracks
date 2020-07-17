function [ mapug, mappg ] = mapping( dmulocal,dmplocal, dmglobal )
%MAPping: a mapping array to set up one-one correspondence between index in
%u array and index in u-p array
%   Detailed explanation goes here
if dmulocal+dmplocal~=dmglobal
    error('The mapping function does not accept such dimension combination')
end
udof=2;
pdof=1;
dof=udof+pdof;
mapug=zeros(1,dmulocal);
% mappg=zeros(1,dmplocal);
for iudof=1:udof
mapug(iudof:udof:end)=iudof:dof:dmglobal;
end
mappg=dof:dof:dmglobal;
end

