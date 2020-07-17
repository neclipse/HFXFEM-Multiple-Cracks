function [ uind, pind  ] = mapping2( upindex )
%MAPping: find uind and pind from the global index in u-p array
%   Detailed explanation goes here
udof=2;
pdof=1;
dof=udof+pdof;
quo=floor(upindex/dof);
rem=mod(upindex,dof);
pflag=(rem==0);
uflag=~pflag;
uind=udof*quo(uflag)+rem(uflag);
pind=quo(pflag);
end

