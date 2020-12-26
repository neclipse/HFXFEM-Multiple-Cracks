function load=ifstd( obj, newmark)
    u=obj.Un2i;      % elemental iterative total displacement vector, u_sup(n+1)_sub(i+1)
    p=obj.Pn2i;      % elemental iterative total pore pressure vector, p_sup(n+1)_sub(i+1)
    v=obj.Un2t1;
    a=obj.Un2t2;
    pd=obj.Pn2t1;
    %% update stress at GaussPntDictM, obsolete after 12/22/20
    % moved to calstress
    a1=newmark.a1;
    % calculate internal loadu, -a1(M*a+K*u-Q*p)
    loadu=-a1*(obj.M*a+obj.K*u-obj.Q*p);
    % calculate internal loadp, Q'*v+H*p+S*pd
    loadp=transpose(obj.Q)*v+obj.H*p+obj.S*pd;
%     loadp=zeros(1,length(p));
    %% need rearrange loadu and loadp to adapt to the data structure[f1x,f1y,f2p]
    load=zeros(length(obj.Locarray),1);
    load(1:3:end-2)=loadu(1:2:end-1);
    load(2:3:end-1)=loadu(2:2:end);
    load(3:3:end)=loadp(:);
%     obj.IntLoadVec=load;            % comprehensive internal load vector,[f1x,f1y,f2p]
% 	obj.IntLoadVec=struct(loadu,loadp);  % element relative internal force vector for current accumulative load increment
end

