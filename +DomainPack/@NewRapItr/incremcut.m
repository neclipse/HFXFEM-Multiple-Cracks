function incremcut( obj )
%INCREMCUT Method of class NewRapItr
%   Cut the current load increment factor into two halves
fprintf('Incremcut is activated at No.%d increment\n',obj.IInc);
iinc=obj.IInc;
if iinc==1
    new=0.5*obj.Increments(iinc);
    obj.Increments=[new,obj.Increments];
elseif iinc>1
    step=0.5*(obj.Increments(iinc)-obj.Increments(iinc-1));
    new=obj.Increments(iinc-1)+step;
    front=obj.Increments(1:iinc-1);
    back=obj.Increments(iinc:end);
    obj.Increments=[front,new,back];
end
end

