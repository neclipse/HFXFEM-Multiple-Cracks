classdef GaussPnt_LE_UP
   properties
       Xi          % s-coordinate in local parent coordinate system
       Eta         % t-coordinate in local parent coordinate system
       H           % integration weight
       X           % x-coordinates in global cartesian coordinate system
       Y           % y-coordinates in global cartesian coordinate system
       NNodes      % number of Nodes involved in the calculation of gaussian points
       GBINP       % Gloabal inputs required 
       Bmat        % B matrix for this gaussian point
       Bmatenr     % a cell array, same structure as the Enf
       DetJ        % Determinant of the Jacobian matrix
       Nu          % Displacement_Shape function for the associated element
       Nuenr       % a cell array, same structure as the Enf
       Npenr       % a cell array, same structure as the Enf
       Np          % pressure_shape function for the associated element
       DNp         % derivative of Np with respect to global coordinate system
       DNpenr      % % a cell array, same structure as the Enf
       Stressp     % current effective stress at this gaussian point
       P           % current pressure at this gaussian point
       Stress      % current total stress at this gaussian point
       Strainc     % relative strain comparing to the initial state
       Uy
       Enf=cell(1,3);         % cell array of struct array{enfid1,enfid2...}, changed on 11/01/2018
%        Tangent_coh % cohesive tangent matrix describing the cohesive crack separation law
%        CrackOpening
%        CrackPerm   % Aparrent Permeability of the crack, depending on the crack opening. 
       % enfid itself is a struct
       % arrary{'Phi',Phi,'Phishift',Phishift,'Ridge',Ridge}
   end
   
   methods
       function obj=GaussPnt_LE_UP(xi,eta,h,nnodes,gbinp)
           if nargin>0
               obj.Xi=xi;
               obj.Eta=eta;
               obj.H=h;
               obj.NNodes=nnodes;             % number of nodes involved in the calculation of gaussian points
               obj.GBINP=gbinp;
           end
       end
       obj=matsu(obj,disp,dp);
       obj=matsu_enriched(obj,us,ue,ps,pe);
       obj=psu(obj);
       obj=preparing(obj,X,Y,EnrichNum);
       obj=enriching(obj);                  % not implemented here though. See Enr.enrichgauss.
       matct(obj);
   end
    
end