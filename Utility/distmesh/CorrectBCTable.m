 function BCType = CorrectBCTable(EToV,BCType,mapnodes,BCcode)
 
 % function BCType = CorrectBCTable(EToV,BCType,mapnodes,BCcode);
 % Purpose: Setup BCType for boundary conditions in 2D
 % By Allan P. Engsig-Karup
 
%  Globals2D;
 
 VNUM = [1 2;2 3;3 1]; % face orientations
 
 for k = 1:K    
     % Test for each edge
    for l = 1:Nfaces 
         m = EToV(k,VNUM(l,1)); n = EToV(k,VNUM(l,2));
 
         % if both points are on the boundary then it is a boundary face!
         ok=sum(ismember([m n],mapnodes));
         if ok==2 
           BCType(l,k)=BCcode; 
         end;
 
    end
end
