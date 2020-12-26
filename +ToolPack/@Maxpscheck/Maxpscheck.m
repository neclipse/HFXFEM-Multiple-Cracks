classdef Maxpscheck < ToolPack.PropagationCheck
% Check the maximum principal stress against the cohesive strength of the
% rock matrix (for the current stage, do not bring up the strength of the
% tensile strings) 06052019
properties
    Stressp     % effective stress at the crack tip
    Elemahead  
    Unstable
    Tol2        % Upper Tolerance for the unstable crack growth
    Dimensionless
    Theta
    Growdirection;
end
methods
    function obj=Maxpscheck(varargin)
        p=inputParser;
        checks= @(x) length(x)==4 && isnumeric(x);
        defaultstressp=zeros(4,1);
        checke= @(x) isa(x,'FEPack.Elem_2d_UP');
        defaultcheckposition='tip';
        validposition={'tip','center','all'};
        checkp= @(x) any(validatestring(x,validposition));
        checkthreshold=@(x) isnumeric(x) && x~=0;
        addRequired(p,'threshold',checkthreshold);
        addOptional(p,'stressp',defaultstressp, checks);
        addOptional(p,'elemahead',[], checke);
        addOptional(p,'checkposition',defaultcheckposition,checkp);
        addParameter(p,'tolerance',0.05);
        addParameter(p,'tolerance2',0.1);
        parse(p,varargin{:});
        obj.Threshold=p.Results.threshold;
        obj.Stressp=p.Results.stressp;
        obj.Elemahead=p.Results.elemahead;
        obj.Mode=p.Results.checkposition;
        obj.Tol=p.Results.tolerance;
        obj.Tol2=p.Results.tolerance2;
        obj.Type='Maxps';
    end
    
    
    % To calculate the principal variable for the crack growth
    function cal_pvariable(obj,stressp,tip,omega,varargin)
        % elem is the element ahead of the crack tip, stressp is the stress
        % calculated at the crack tip EnrCrackTip has access to this info.
        % Call obj.GrowCheck.cal_pvariable(obj.elementahead,obj.Stressp)
        if ~isempty(varargin)
            obj.Elemahead=varargin{1};
        end
        obj.PVariable=[];
        obj.Theta=[];
        obj.Stressp=stressp;    % nonlocal effective stress at the tip
        obj.psu_exact(obj.Stressp,tip,omega);
        % Note the stress values at guassian points have been
        % calculated in the call of tip.calstress_nonlocal.
        switch obj.Mode
            case 'tip'
                % do nothing 
            case 'center' % not recommended for enriched elemahead.
                % For current stage, the Elemahead is always standard element.
                % Update: the elemahead can be an enriched element.12/22/20
                if obj.Elemahead.EnrichNum==0
                    % the stressp at the element center is already calculated
                    % in the tip.calstress_nonlocal.
                    obj.psu_exact(obj.Elemahead.Stressp,tip,omega);
                end
                % % Need to think about if center stress is appropriate for
                % element containing a crack.
            case 'all' % not recommended for enriched elemahead.
                % check the gaussian points within the element ahead of the
                % tip 
                if obj.Elemahead.EnrichNum>0
                    for ig=1:length(obj.Elemahead.EnrichGauss)
                        gaussstressp=obj.Elemahead.EnrichGauss(ig).Stressp;
                        obj.psu_exact(gaussstressp,tip,omega);
                    end
                else
                    for ig=1:length(obj.Elemahead.GaussPntDictM)
                        gaussstressp=obj.Elemahead.GaussPntDictM(ig).Stressp;
                        obj.psu_exact(gaussstressp,tip,omega);
                    end
                end

        end
    end
    % Check the Pvariable against the threshold (cohesive strength of matrix) 
    function growcheck(obj)
        % The growcheck need to be fine tuned based on further observations
        % 06172019, because it is yet clear how the modes will affect the
        % check.
        %         obj.cal_pvariable; % called in EnrCrackTip separately
        f=obj.PVariable/obj.Threshold;
        obj.Dimensionless=f;
        % note the check should be highly dependent on the obj.Mode. 
        % UPDATED ON 08132019 WHEN IT IS NOTICED THAT THESE CRITERIA HIGHLY
        % AFFECT THE CRACK PROPAGATION.
        switch obj.Mode
            case 'tip'
                if f>1+obj.Tol2
                    %             warning('cut back the current time increment');
                    obj.Unstable=true;
                    obj.Growflag=false;
                elseif f>=1+obj.Tol
                    %             warning('cut the following increments');
                    obj.Unstable=true;
                    obj.Growflag=true;
                elseif f>1
                    obj.Unstable=false;
                    obj.Growflag=true;
                else
                    obj.Unstable=false;
                    obj.Growflag=false;
                end
            case 'center'
                if max(f)>1+obj.Tol2
                    %             warning('cut back the current time increment');
                    obj.Unstable=true;
                    obj.Growflag=false;
                elseif max(f)>=1+obj.Tol 
%                     warning('cut the following increments');
                    obj.Unstable=true;
                    obj.Growflag=true;
                elseif max(f)>1
                    % Here,as long as the pvariable at the center pass the
                    % threshold, the criterion is met.
                    obj.Unstable=false;
                    obj.Growflag=true;
                else
                    obj.Unstable=false;
                    obj.Growflag=false;
                end
            case 'all'
                temp=f(f>1);
                if max(f)>1+obj.Tol2
                    %             warning('cut back the current time increment');
                    obj.Unstable=true;
                    obj.Growflag=false;
                elseif max(f)>=1+obj.Tol && sum(f)/length(f)>=1
                    %             warning('cut the following increments');
                    obj.Unstable=true;
                    obj.Growflag=true;
                    %             obj.Growdirection=sum(obj.Theta)/length(obj.Theta);
                elseif length(temp)>2
                    obj.Unstable=false;
                    obj.Growflag=true;
                else
                    obj.Unstable=false;
                    obj.Growflag=false;
                end
        end
    end
    function theta = growdirection(obj)
        % theta should be determined based on the obj.mode and f value 
        theta=obj.Theta(1); % Based on the nonlocal tip stress
%         switch obj.Mode
%             case 'tip'
%                 % check the tip only
%                 theta=obj.Theta(1);
%             case 'center'
%                 temp=obj.Theta(obj.Dimensionless>1);
%                 theta=sum(temp)/length(temp);
%             case 'all'
%                 % check the gaussian points within the element ahead of the tip
%                 temp=obj.Theta(obj.Dimensionless>1);
%                 theta=sum(temp)/length(temp);
%         end
        % theta is the angle between the current global axes to the
        % principal global axes if positive then counterclockwisely
        % rotated if negative then clockwisely rotated
        %IMPORTANT BUG: DO NOT USE THETA MIMUS OBJ.OMEGA. 08082019
        % ALSO USE ABS(THETA)
        if abs(theta)<1e-5
            theta=0;
        end
        obj.Growdirection=theta;
%       11/10/2019 for simplicity to debug. 
%         theta=0;
        obj.Growdirection=theta;
    end
    [psmax,gtheta] = psu_exact(  obj, s,  tip, omega );

end   
end