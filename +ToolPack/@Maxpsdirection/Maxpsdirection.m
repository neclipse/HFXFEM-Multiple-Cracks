classdef Maxpsdirection < ToolPack.CrackGrowDirectionLaw
% Chang Huang, the most basic version for grow direciton. 06062019
   methods
       function obj=Maxpsdirection(omega)  
           obj.Omega=omega;
            % The effective stress at the tip is passed in as the PVariable
       end
       
       function theta = growdirection(obj,stressp)
            obj.PVariable=stressp;
            [~,theta]=psu(obj.PVariable);
            % theta is the angle between the current global axes to the
            % principal global axes if positive then counterclockwisely
            % rotated if negative then clockwisely rotated
            %IMPORTANT BUG: DO NOT USE THETA MIMUS OBJ.OMEGA. 08082019
            % ALSO USE ABS(THETA)
            % changed back on 0816: should use theta directly.
%             if abs(theta)<1e-6
%                 theta=0;
%             end
            obj.Theta=theta;
       end
       
   end
    
    
end