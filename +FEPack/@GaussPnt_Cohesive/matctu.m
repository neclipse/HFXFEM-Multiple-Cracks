function  [traction,stagechangeflag,obj]=matctu(obj,ua,Due)
%MATCTU cohesion traction Update method of GaussPnt_Cohesive (value calss)
% called by elem_2d_UP.ifstd_enriched
% update crackdisp within obj, note this obj is not returned to the obj
% so that the original crack opening and other information is not updated
% from here
% obj.updatecrackopening(us,ue);

%% 07302019, A fully Newton_Raphson scheme is required to handle the
% complicated TSL. The current form is only okay as of monotonic loading or
% unloading. As long as the cohesive tangent matrix changes abruptly from
% during calcultion, it is hard to reach convergence. Here, no matter how
% the cohesion is updated, the Kc in the linear equaiton system is not
% changed. That is the problem of the modified Newton-Raphson scheme.
% However, it is usually the case of monotonic loading for hydraulic
% fracturing process so the scheme is not changed to full NR scheme yet.

%%
if obj.InitialMode>2 % matctu is not needed for perforated or smeared crack
    % ul=obj.Amat*obj.CrackDisp;     % ul is local displcaement discontinuity averaged from nodal values
    %% Part 2:Derive Tangent_loc in the local orthogonal coordinate system
    % obj.TractionLaw.Disp=ul;            % Local
    % if isempty(varargin)
    % obj.TractionLaw.checkseparation(2); % check trial separation stage.
    stagechangeflag=false;
    % reset the AfterPeak at every iteration.
    obj.TractionLaw.AfterPeak=false;
    if obj.TractionLaw.Separated==1
        traction=[0;0];
        return;
    end
    traction=obj.TractionO+obj.Tangent_coh*(obj.Nuenrplus-obj.Nuenrminus)*Due;
    %% Validate or modify the traction
    % crack opening in the last increment
    %     localdisp=obj.Amat*obj.CrackDisp;
    %     lastopening=localdisp(2);
    %     initialdisp=obj.Amat*obj.IniCrackDisp;
    %     initialopening=initialdisp(2);
    % calculate the current crack opening
    % the two ways: (Nuenrplus-Nuenrminus)*ua or Nu*ua
    % choose one and keep consistent in this function. 10/22/20
    CrackDisp=(obj.Nuenrplus-obj.Nuenrminus)*ua+obj.IniCrackDisp;
    % 12/20/20 What if I update crack disp here? Does it make matct in
    % later iterations calculated based on the crackdisp from this
    % iteration? No, it does not. because matctu does not return obj.
    obj.CrackDisp=CrackDisp; 
    calcrackopening=obj.Ntaud'*CrackDisp;
    %     shearopening=obj.Mtaud'*CrackDisp;
    %% By current opening, change the stiffness matrix also
    if calcrackopening<=0 
        Tangent_lol=[obj.TractionLaw.Tshearc,0;0,obj.TractionLaw.Tnormalc];
        Tangent_coh=obj.Amat'*Tangent_lol*obj.Amat;
        if any(ua)  % to avoid disturbing the initial call of intforcer when ua are all zeros.
            temp1=Tangent_coh*(obj.Nuenrplus-obj.Nuenrminus)*ua;
            temp2=obj.Amat*temp1; % from global x-y to local shear-normal
            tnormal=temp2(2);
            % Bug: traction(1) is global traction in x-direction, it may
            % not be tshear. Change it to temp2(1)
            tshear=temp2(1);
            if tnormal>0  % normal traction becomes tensile
                tnormal=obj.TractionLaw.IniTraction-tnormal;
                traction=[0;tnormal];
                traction_eff=sqrt(tshear^2+tnormal^2);
                if traction_eff>obj.TractionLaw.PeakTraction*1.1
                    %                     release the shear traction if tnormal turns positive
                    traction=obj.Amat'*[0;obj.TractionLaw.IniTraction];
                end
            else
                % compressive mode, could use mohr-coulomb criterion, not
                % implemented yet. 12/19/20.
                traction=temp1;
            end
        end
    else % current opening is positive
        % if the previous increment is already stage 1
        if obj.TractionLaw.Stage==1 % only need to worry about stage change at stage 1
            traction_eff=sqrt(traction(1)^2+traction(2)^2);
            traction_effold=sqrt(obj.TractionO(1)^2+obj.TractionO(2)^2);
            traction_inc=abs(traction_eff-obj.TractionLaw.IniTraction);
            traction_incold=abs(traction_effold-obj.TractionLaw.IniTraction);
            legitmate_inc=abs(obj.TractionLaw.KrgTraction-obj.TractionLaw.IniTraction);
            if traction_inc>legitmate_inc*1.15 && legitmate_inc>0.2*obj.TractionLaw.PeakTraction
                % For combined traction law
                if  traction_incold<legitmate_inc*0.85
                    % this 0.99 ratio is to avoid cutting the increment too many
                    % times when the TractionO is around the peaktraction.
                    stagechangeflag=true;
                    return;
                else
                    obj.TractionLaw.AfterPeak=true;
                    % stagechangeflag=true;
                    traction=obj.Amat'*obj.TractionO;   % Just use the old traction
                    return;
                end
            end
            % if the previous increment is still at stage -1
        elseif obj.TractionLaw.Stage==-1
            traction_eff=sqrt(traction(1)^2+traction(2)^2);
            if traction_eff>obj.TractionLaw.PeakTraction
                stagechangeflag=true;
                return;
            end
        end
    end
else
    stagechangeflag=false;
    traction=[0;0];
end
% Do not assign tractio to obj.Traction here as it can be temporarty for
% stagecheck.
end

