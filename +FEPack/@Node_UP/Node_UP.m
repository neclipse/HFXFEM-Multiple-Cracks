classdef Node_UP < handle
   properties
       Ind
       Type
       X
       Y
       DofArray      % global dof index array
       EnrDofArray=cell(1,10);  % global enrdof index array
       UArray        % global disp dof index array
       UenrArray=cell(1,10);     % a cell array, same structure as the Enf
       PArray        % global pressure dof index array
       PenrArray=cell(1,10);     % a cell array, same structure as the Enf
       pressure      % calculated pressure by the linear system creater
       ItrDispList   % Iterative values of dofs in corresponding order as DofArray, column array
       IncDispList   % Incremental values of dofs in corresponding order as DofArray,column array
       TotDispList   % Accumulated values of dofs in corresponding order as DofArray,column array
       NoDofs
       NoUDofs
       NoPDofs
       NoEnrDofs=zeros(1,10);    %  array, same structure as the Enrich
       NoUenrDofs=zeros(1,10);    
       NoPenrDofs=zeros(1,10);
       Stress=zeros(5,1)         % Total stresses vector [repetition,sx,sy,sxy,sz]
       Stressp=zeros(5,1)        % Effective stresses vector [rep,sx,sy,sxy,sz]
       Leakoff=0;                % Leakoff rate
       ALeakoff=0;               % Accumulated leakoff volume
       CInjection=0;             % Concentrated injection to this node, calculated from encrack.cal_qextenr
       Enrich=zeros(1,10);       % Enrich item id,[nan,id1,id2,...]
   end
   
   methods
       function obj= Node_UP(ind,type,x,y)
           if nargin>0
           obj.Ind=ind;
           obj.Type=type;
           obj.X=x;
           obj.Y=y;
           obj.NoUDofs=2;
           obj.NoPDofs=1;
           obj.NoDofs=3;
           end
       end
%        function NoPDofs=get.NoPDofs(obj)
%            if obj.Type==1
%                NoPDofs=1;
%            elseif obj.Type==2
%                NoPDofs=0;
%            end
%        end
%        function NoDofs=get.NoDofs(obj)
%                NoDofs=obj.NoUDof+obj.NoPDofs;
%        end
        function setenrich(obj,id)
            % obj.Enrich stores all the ids of involved enrichitem with
            % this element
            L=id~=obj.Enrich;
            if all(L)
                obj.Enrich(id)=id;
            end
        end
       function [totndof,totudof,totpdof]=givedofarray(obj,totndof,totudof,totpdof)
           obj.DofArray=totndof+1:totndof+obj.NoDofs; % row array
           obj.UArray=totudof+1:totudof+obj.NoUDofs;
           obj.PArray=totpdof+1:totpdof+obj.NoPDofs;
           totndof=totndof+obj.NoDofs;
           totudof=totudof+obj.NoUDofs;
           totpdof=totpdof+obj.NoPDofs;
       end

       % enriched dof array is assigned by the Enf.addnodedofs and EnrichItem.updatedofarray function
       % because different ENF may assign different number of dofs and
       % there may exist multiple EnrichItems
%        function [totnenrdof, totuenrdof, totpenrdof]=giveenrdofarray(obj,totnenrdof, totuenrdof, totpenrdof)
%            % separate standard dofs and enriched dofs in global indexing
%            
%        end
   end
end