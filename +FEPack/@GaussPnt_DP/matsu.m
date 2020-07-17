function  matsu(obj)
%MATSU Stress State Update method of GaussPnt_DP
%   Elastic predictor/ Plastic returning map method is employed
%   For every iteration, this method is called to calculate the accumulated stress for current load increment
%   TENSILE STRESS IS DEFINED POSITIVE
%% Preparing       
% retrieve material properties from element class
global  G K xi eta etabar c0 H Delastic inistress;             % initial stress state
% initialize some parameters for this return
epflag=0;
apexflag=0;
epbar=obj.EPBar;                            % epbar for current iteration, epbar_sub{n+1}
tds=zeros(4,1);                             % Elastic trial deviatoric stress
stress=zeros(4,1);                          % Updated current total stress
tol=1e-4;                                   % returning tolerance
%% Elastic predictor
% Last returned total elastic strain tensor + strain change at this iteration
etrial=obj.ElaStn+obj.ItrStn;              
etv=etrial(1)+etrial(2)+etrial(4);          % Relative elastic trial volumetric strain;
tps=3*K*etv/3;                                  % Trial Pressure stress
% -Trial deviatoric stress
tds([1,2,4])=2*G*(etrial([1,2,4])-etv/3);
tds(3)=G*etrial(3);                         % G*gama=tau; Elastic trial strain is in engineering strain form. 
% Compute elastic trial stress J2 invariant
SQRJ2T=sqrt(tds(3)^2+0.5*(tds(1)^2+tds(2)^2+tds(4)^2));
cohe=c0+H*epbar;
	%% Check plastic admissity
	phi=SQRJ2T+eta*tps-xi*cohe;             % Drucker-Prager yield surface function
    if cohe~=0
		residual=phi/abs(cohe);
	end
	if residual>tol
		epflag=1;
	end
	%% Plastic Return mapping
    if epflag==1 
        %Return to the smooth cone first
      % Solve the plastic multiplier in closed form for linear
      % hardening Drucker-Prager model. The formula also applies to perfect
      % plastic model, for which H equals to zero.
        dgama=phi/(G+K*eta*etabar+H*xi^2); 
        SQRJ2=SQRJ2T-G*dgama;   % Updated sqrt(J2_sub{n+1}), not trial value
        if SQRJ2>0
            factor=1-G*dgama/SQRJ2T;
        elseif SQRJ2==0
            factor=0;        
        elseif SQRJ2<0
            apexflag=1;
        end
        % return to the apex
        if apexflag==1
            depv=(tps-cohe*xi/etabar)/(xi^2*H/(eta*etabar)+K);
            dgama=depv/etabar;
            factor=0;
        end
        %% Update stress, strain and epbar
        % -Update total Stress
        %% WRONG UPDATED TOTAL STRESS
        ps=tps-K*dgama*etabar;
        stress([1,2,4])=factor*tds([1,2,4])+ps;             
        stress(3)=factor*tds(3);
        stressc=stress-inistress;
        obj.Stressc=stressc;
        % -Update epbar and dgama
        obj.EPBar=epbar+xi*dgama;
        obj.DGama=dgama;
        % -Update elastic (engineering) strain components
        elastn=Delastic\stress;
        delastn=elastn-obj.ElaStn;      % Change in elastic strain
        obj.ElaStn=elastn;              % Updated elastic strain
        % Save relative elastic trial strain for the calculation of consistent tangent
        % operator
        obj.ETrial=etrial;
        % Update relative plastic strain components
        dplastn=obj.ItrStn-delastn;     % Change in plastic strain
        obj.PlaStn=obj.PlaStn+dplastn;  % Updated plastic strain
        %% Pure Elastic step
        elseif epflag==0
        stress([1,2,4])=tds([1,2,4])+tps;
        stress(3)=tds(3);
        stressc=stress-inistress;
        obj.Stressc=stressc;
        % Update relative strain
        obj.ElaStn=etrial;
        obj.ETrial=etrial;
    end
% Update some algorithmic variables before exit
obj.EPFlag=epflag;
obj.ApexFlag=apexflag;
end

 