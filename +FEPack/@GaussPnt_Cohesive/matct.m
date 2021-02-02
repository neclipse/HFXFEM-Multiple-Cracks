function  obj=matct(obj,varargin)
%MATCT Method of GaussPnt_Cohesive
%Calculate the Cohesive tangent matrix in local coordinate system and
% then transform it to global coordinate system
% Note ua is the displacement discontinuity in global coordinate system,
% the cohesive tangent matrix is a function of the ua according to the
% cohesive crack model.
%% Part 1:Transformation matrix A in Khoei 2012,from the global to the local
if obj.InitialMode>2 % matct is not needed for perforated or smeared crack
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
    ul=Amat*obj.CrackDisp;     % ul is local displcaement discontinuity averaged from nodal values
    %% Part 2:Derive Tangent_loc in the local orthogonal coordinate system
    obj.TractionLaw.Disp=ul;            % Local [dshear, dnormal]
    % if isempty(varargin)
    % obj.TractionLaw.checkseparation (1); % trial separation stage.
    if obj.TractionLaw.AfterPeak
        Tangent_loc=obj.TractionLaw.givepeaktangent;
    else
        Tangent_loc=obj.TractionLaw.givetangent;
    end
    %% Part 3:Transform Tangent_loc form local to global, T_global=lamda'*T_loc*lamda
    obj.Tangent_coh=Amat'*Tangent_loc*Amat;
    % else
    %     stage=obj.TractionLaw.checkseparation;
else
    obj.Tangent_coh=[0,0;0,0];
end

