function crtRHS(obj)
%CRTRHS Create the right hand side of the linear system
%   This function include integration and assembly
FFG=zeros(obj.Drow,1);
nonodes=obj.ElemDict(1).NoNodes; %   nonodes: the number of nodes invloved in the element stiffness matrix
[~,Nline_xi]=lineshape(nonodes,1);
n=length(Nline_xi); % An indirect way to decide the number of nodes along a side of element
noloads=size(obj.Loads,1); % total number of boundaries
if nonodes==3 || nonodes==6 % For triangular elements
    [xx,ww]=GaussQuad(1,1);% calculating gauss points along a line from 0 to 1
elseif nonodes==4 || nonodes==8 % For rectangular elements
    [xx,ww]=GaussQuad(1);% calculating gauss points along a line from -1 to 1
end
% xx is the coordinate, and ww is the correponding weight 
%% Calculate element force vetor fe boundary by boudary
ind1=logical(obj.BCTablel(:,end)==2); % find all rows where the boundary condition for loads equals 2
BCTable=obj.BCTablel(ind1,:);
% loop over boundaries
for inb=1:noloads
   ind=find(BCTable(:,1)==inb); % find all the rows where the first column equals inb
   table=BCTable(ind,:); % Extract one boundary condition table from the comprehensive boundary condition table
   load=obj.Loads(inb,:)';
   % loop over one boundary
   for i=1:length(ind)
      fe=zeros(2*n,1);
      locarray=zeros(2*n,1);
      X=[obj.NodeDict(table(i,2:(1+n))).X]; % X is a row vector
      % X is the x coordinates list of the nodes along the boundary side
      Y=[obj.NodeDict(table(i,2:(1+n))).Y]; 
      % loop over gaussian points
      for j=1:length(xx)
         xi=xx(j);
         [Nline,Nline_xi]=lineshape(nonodes,xi);
         x_xi=Nline_xi*X';
         y_xi=Nline_xi*Y';
         % THIS FORMULA IS DEBATABLE for Triangle Elements, because it is from [0,1] which should be changed [-1,1] (02152019)
         normal=[y_xi,0;0,-x_xi];    
         fet=Nline'*normal*load;     % THIS FORMULA IS DEBATABLE for Triangle Elements
         fe=fe+ww(j)*fet;
      end
%% Assembly into global force vector
      for k=1:n
         locarray(2*k-1:2*k)=obj.NodeDict(table(i,1+k)).DofArray(1:2);
         % using doflist in the node object to construct a locarray for
         % assembly
      end
      FFG(locarray)=FFG(locarray)+fe;
   end
end
obj.RHS=FFG;
end

