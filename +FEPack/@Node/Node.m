classdef Node < handle
   properties
       Ind
       Type
       X
       Y
       DofArray      % global dof index array
       pressure      % calculated pressure by the linear system creater
       ItrDispList   % Iterative values of dofs in corresponding order as DofArray, column array
       IncDispList   % Incremental values of dofs in corresponding order as DofArray,column array
       TotDispList   % Accumulated values of dofs in corresponding order as DofArray,column array
       
   end
   properties(Dependent)
       NoDofs
   end 
   methods
       function obj= Node(ind,type,x,y)
           if nargin>0
           obj.Ind=ind;
           obj.Type=type;
           obj.X=x;
           obj.Y=y;
           end
       end
       function NoDofs=get.NoDofs(obj)
          if obj.Type==1    
              NoDofs=2;
          elseif obj.Type==2
              NoDofs=2;
          end
       end
       function totndof=givedofarray(obj,totndof)
           obj.DofArray=zeros(1,obj.NoDofs); % row array
           for idof=1:obj.NoDofs
           totndof=totndof+1;  
           obj.DofArray(idof)=totndof;
           end
       end
   end
end