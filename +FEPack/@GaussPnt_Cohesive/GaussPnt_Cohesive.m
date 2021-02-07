classdef GaussPnt_Cohesive < FEPack.GaussPnt_LE_UP
   properties
       Tangent_coh      % cohesive tangent matrix describing the cohesive crack separation law
       CrackDisp;       % crack displacement at current gauss point in global coordinate system
       IniCrackDisp;
       MinCrackOpening  % for perforated initial notch, 10282019 also used to compare with Abaqus (2e-3 m)
       CrackOpening;
       Uplus            % displacement at the positive face of the crack, positive means the signed distance function is positive, outdated 10/20/20,
       Uminus           % displacement at the negative face of the crack, they are rigorous again as full Nuenr is recovered due to new chagnes 11/27/20 
       Nuenrplus        % enriched shape function right above the crack Nuenr+
       Nuenrminus       % enriched shape function right beneath the crack Nuenr-
       FractureP        % Fracture pressure interpolated from standard pdof and all penrdofs at the element nodes
       Ds               % The interval length for integral
       Ntaud            % unit normal vector
       Mtaud            % unit tangent vector
       Amat             % 2d-Coordinate transformation matrix from global to the local;
       TractionO;       % Converged traction 
       Traction;        % Traction (MPa), [tx,ty] 
       TractionLaw      % An object of TractionLaw
       InitialMode      % Same flag as inherited from EnrCrackBody, but may change for newly created tensile fractures.
       Alpha=pi/2;            % the angle between the loading and the crack plane for
                        % initially inplace mode, counterclockwise from plane to
                        % loading.
   end
   
   methods
       function obj=GaussPnt_Cohesive(xi,eta,h,nnodes,gbinp,initialmode,alpha)
           if nargin==0
               super_args={};
           elseif nargin>=5
               super_args={xi,eta,h,nnodes,gbinp};
           else
               error('Wrong number of inputs')
           end
           obj=obj@FEPack.GaussPnt_LE_UP(super_args{:});
           if nargin>=5
               obj.InitialMode=initialmode;
               obj.Alpha=alpha;
           end
       end
       
       function obj=initiate(obj,varargin) 
          % initialize the crack displacement and traction 04042019, was
          % not adopted in the latest version 06042019.
          if ~isempty(varargin)
              minaperture=varargin{1};
              perfapeture=varargin{2};
          end
          Vx=[1,0];
          Vy=[0,1];
          Vxl=obj.Mtaud;
          Vyl=obj.Ntaud;
          % Refer to the coordinate transformation
          lxl=Vx*Vxl;     % cos(theta)
          mxl=Vy*Vxl;     % sin (theta)
          lyl=Vx*Vyl;     % -sin(theta)
          myl=Vy*Vyl;     % cos(theta)
          % 2d-Coordinate transformation matrix from global to the local;
          obj.Amat=[lxl,mxl;lyl,myl];
          switch obj.InitialMode
              case 1 % perforated mode
                  obj.MinCrackOpening=perfapeture;
                  obj.CrackOpening=obj.MinCrackOpening;
                  obj.CrackDisp=obj.Amat'*[0;0];
                  obj.IniCrackDisp=obj.CrackDisp;
                  obj.Traction=[0;0];
                  obj.TractionO=obj.Traction;
              case 2 % smeared mode
                  obj.MinCrackOpening=0;
                  obj.CrackOpening=0;
                  obj.CrackDisp=obj.Amat'*[0;0]; % initially 
                  obj.IniCrackDisp=obj.CrackDisp;
                  obj.Traction=[0;0];
                  obj.TractionO=obj.Traction;
              case 3 % compressive mode
                  initraction=0; % will be updated after first increment
                  ul=[0;-1e-16]; % [us;un], local displcaement discontinuity averaged from nodal values
                  obj.CrackDisp=obj.Amat'*ul;
                  obj.IniCrackDisp=obj.CrackDisp;
                  obj.MinCrackOpening=0;    % Abaqus setttings
                  obj.CrackOpening=obj.Ntaud'*obj.CrackDisp; % should be un
                  obj.Traction=obj.Amat'*[initraction*cos(obj.Alpha);initraction*sin(obj.Alpha)];
                  obj.TractionO=obj.Traction;
              case 4 % tensile mode for newly propagated crack segment
                  % The initraction can also be directly calculated from the
                  % stress states at the linegauss points given the linegauss
                  % points are predefined. 10/06/20 if developed from smeared crack
                  % THE SWITCH BLOCK NEED TO BE UPDATED AFTER THE SMEARED CRACK 01/31/2021
                  switch obj.TractionLaw.Type
                      case 'linear'
                          initraction=(1-obj.TractionLaw.Lambdaini)*obj.TractionLaw.PeakTraction;
                      case 'unified'
                          initraction=obj.TractionLaw.IniTraction+obj.TractionLaw.Lambdaini*...
                              (obj.TractionLaw.PeakTraction-obj.TractionLaw.IniTraction)/obj.TractionLaw.Lambdacr;
                      case 'bilinear'
                          initraction=obj.TractionLaw.Lambdaini/obj.TractionLaw.Lambdacr*obj.TractionLaw.PeakTraction;
                  end
                  separation=obj.TractionLaw.Lambdaini*obj.TractionLaw.CriticalDisp;
                  ul=[separation*cos(obj.Alpha);separation*sin(obj.Alpha)]; % [us,un], local displcaement discontinuity averaged from nodal values
                  obj.CrackDisp=obj.Amat'*ul;
                  obj.IniCrackDisp=obj.CrackDisp;
                  obj.MinCrackOpening=minaperture;    % Abaqus setttings
                  obj.CrackOpening=separation*sin(obj.Alpha)+obj.MinCrackOpening;
                  obj.Traction=obj.Amat'*[initraction*cos(obj.Alpha);initraction*sin(obj.Alpha)];
                  obj.TractionO=obj.Traction;
                  % BUG:SHOULD NOT CHANGE THE TRACTIONLAW.INITRACTION AS IT
                  % IS USED TO CALCULATE THE TANGENT. 02/01/2021
%                   obj.TractionLaw.IniTraction=initraction;
          end

       end
       function obj=updatecrackopening(obj,us,ua,varargin) 
           % outdated Uplus and Uminus 10/20/20 because they are not highly
           % necessary and only used in domain.snapshot and
           % crackbody.postprocess. The important thing is the
           % obj.CrackDisp and obj.CrackOpening. The simplified Nuenrplus
           % and Nuenrminus are good for these calculation but not Uplus
           % and Uminus. 
           % reverted 11/27/20, Uplus and Uminus are not outdated.
           if obj.InitialMode~=2
               obj.Uplus=obj.Nu*us+obj.Nuenrplus*ua;
               obj.Uminus=obj.Nu*us+obj.Nuenrminus*ua;
               % add the obj.IniCrackDisp as the resultant crackdisp should be
               % used for givetangent.
               obj.CrackDisp=(obj.Nuenrplus-obj.Nuenrminus)*ua+obj.IniCrackDisp;
               % it is equivalent to obj.Nu*ua+obj.IniCrackDisp;
               crackopening=obj.Ntaud'*obj.CrackDisp;
               obj.CrackOpening=crackopening+heaviside(crackopening)*obj.MinCrackOpening;
           end
           % one place to prevent interpenetration, 07052019
           % Other places to prevent interpenetration, cohesive and contact
           % behavior. (partially done)
           % IMPORTANT BUG: DO NOT ADD CALCRACKOPENING AGAIN AS
           % OBJ.INICRACKDISP IS ALREADY ADDED TO THE OBJ.CRACKDISP. 091319
       end
       
       %% Prototype
       obj=matsu(obj,us,ps);
       obj=matsu_enriched(obj,us,ue,ps,pe);
       [traction,stagechangeflag,obj]=matctu(obj,ue,Due);
       obj=preparing(obj,X,Y,EnrichNum);       
       obj=enriching(obj);          % preparing the enrichment matrices, Bmatenr,Nuenr,Npenr
       obj=matct(obj);               % calculate the tangent_cohesive in the global coordinate system
   end
    
end