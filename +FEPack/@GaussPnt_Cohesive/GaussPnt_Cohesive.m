classdef GaussPnt_Cohesive < FEPack.GaussPnt_LE_UP
   properties
       Tangent_coh      % cohesive tangent matrix describing the cohesive crack separation law
       CrackDisp;       % crack displacement at current gauss point in global coordinate system
       IniCrackDisp;
       MinCrackOpening  % for perforated initial notch, 10282019 also used to compare with Abaqus (2e-3 m)
       CrackOpening;
       Uplus            % displacement at the positive face of the crack, positive means the signed distance function is positive, outdated 10/20/20
       Uminus           % displacement at the negative face of the crack, outdated, as Nuenrplus and Nuenrminus are not rigorous (only contain partial Nuenr)
       Nuenrplus        % enriched shape function right above the crack Nuenr+, take the same shape of Nuenr, but only contains effective info from the current crack
       Nuenrminus       % enriched shape function right beneath the crack Nuenr-, take the same shape of Nuenr, not very useful as only crackopening is required.
       FractureP        % Fracture pressure interpolated from standard pdof and all penrdofs at the element nodes
       Ds               % The interval length for integral
       Ntaud            % unit normal vector
       Mtaud            % unit tangent vector
       TractionO;       % Converged traction 
       Traction;        % Traction (MPa), [tx,ty] 
       TractionLaw      % An object of TractionLaw
       Perforated       % Initially fully debonded (perforated) 07042019
%        Inplace=false;   % If the gauss_cohesive point belongs to a initially in-place crack
       Alpha      %% the angle between the loading and the crack plane for
                        % initially inplace mode, counterclockwise from plane to
                        % loading.
   end
   
   methods
       function obj=GaussPnt_Cohesive(xi,eta,h,nnodes,gbinp,perforated,varargin)
           if nargin==0
               super_args={};
           elseif nargin>=5
               super_args={xi,eta,h,nnodes,gbinp};
           else
               error('Wrong number of inputs')
           end
           obj=obj@FEPack.GaussPnt_LE_UP(super_args{:});
           if nargin==6
               obj.Perforated=perforated;
           elseif nargin==7
               obj.Perforated=perforated;
               obj.Alpha=varargin{1};
%                if varargin{1}~=pi/2
%                    obj.Inplace=true;
%                end
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
          myl=Vy*Vyl;     % cos(theat)
          % 2d-Coordinate transformation matrix from global to the local;
          Amat=[lxl,mxl;lyl,myl];
          if ~obj.Perforated        % Initially Bonded mode
              % The initraction can also be directly calculated from the 
              % stress states at the linegauss points given the linegauss
              % points are predefined. 10/06/20
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
              ul=[separation*cos(obj.Alpha);separation*sin(obj.Alpha)]; % [us,un]
              obj.CrackDisp=Amat'*ul;     % ul is local displcaement discontinuity averaged from nodal values
              obj.IniCrackDisp=obj.CrackDisp;
              obj.MinCrackOpening=minaperture;    % Abaqus setttings
              obj.CrackOpening=separation*sin(obj.Alpha)+obj.MinCrackOpening;
              obj.Traction=Amat'*[initraction*cos(obj.Alpha);initraction*sin(obj.Alpha)];
              obj.TractionO=obj.Traction;
              obj.TractionLaw.IniTraction=initraction;
          else    % Perforated mode
              obj.MinCrackOpening=perfapeture;
%               obj.CrackOpening=obj.TractionLaw.CriticalDisp;
              obj.CrackOpening=obj.MinCrackOpening;
              obj.CrackDisp=Amat'*[0;0];
              obj.IniCrackDisp=obj.CrackDisp;
              obj.Traction=[0;0];
              obj.TractionO=obj.Traction;
          end

       end
       function obj=updatecrackopening(obj,us,ua,varargin) 
           % update crack opening using the updated crack displacement,
           % not used for the first increment
           % ua is the vector of displacement discontinuity of the all
           % nodes within the current element, zero for standard nodes.
%            if isempty(varargin)
%                % Type 1: standard heaviside
%                obj.CrackDisp=obj.Nu*ua;
%            elseif varargin{1}==1;
%                obj.CrackDisp=obj.Nu*ua;
%            elseif varargin{1}==2;
%                % Type 2: sign-step function
%                obj.CrackDisp=obj.Nu*ua*2;
%            end

           % outdated Uplus and Uminus 10/20/20 because they are not highly
           % necessary and only used in domain.snapshot and
           % crackbody.postprocess. The important thing is the
           % obj.CrackDisp and obj.CrackOpening. The simplified Nuenrplus
           % and Nuenrminus are good for these calculation but not Uplus
           % and Uminus.
           obj.Uplus=obj.Nu*us+obj.Nuenrplus*ua;
           obj.Uminus=obj.Nu*us+obj.Nuenrminus*ua;
           % add the obj.IniCrackDisp as the resultant crackdisp should be
           % used for givetangent.
           obj.CrackDisp=(obj.Nuenrplus-obj.Nuenrminus)*ua+obj.IniCrackDisp;
           % it is equivalent to obj.Nu*ua+obj.IniCrackDisp;
%            obj.CrackDisp=(obj.Nuenrplus-obj.Nuenrminus)*ua;
           obj.CrackOpening=obj.Ntaud'*obj.CrackDisp+obj.MinCrackOpening;
           % one place to prevent interpenetration, 07052019
           % Other places to prevent interpenetration, cohesive and contact
           % behavior. (partially done)
           % IMPORTANT BUG: DO NOT ADD CALCRACKOPENING AGAIN AS
           % OBJ.INICRACKDISP IS ALREADY ADDED TO THE OBJ.CRACKDISP. 091319
%            if calcrackopening>0
%                obj.CrackOpening=obj.MinCrackOpening+calcrackopening;
%            end
       end
       
       %% Prototype
       obj=matsu_enriched(obj,us,ue,ps,pe)
       [traction,stagechangeflag]=matctu(obj,ue,Due);
       obj=preparing(obj,X,Y,EnrichNum);       
       obj=enriching(obj);          % preparing the enrichment matrices, Bmatenr,Nuenr,Npenr
       obj=matct(obj);               % calculate the tangent_cohesive in the global coordinate system
   end
    
end