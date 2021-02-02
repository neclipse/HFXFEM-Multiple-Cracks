function  iterating( obj, iinc, stdpdofs, varargin)
%Iterating Method of NewtonRaphson class
% Loop over equilibrium iterations for one increment, modified
% Newton-Raphson method as the tangent matrix is updated per increment.
% Assemble, solve Linear equation system and update
inp=inputParser;
addOptional(inp,'psddofs',[]);
addOptional(inp,'inclist',[]);
parse(inp,varargin{:});
psddofs=inp.Results.psddofs;
inclist=inp.Results.inclist;
obj.IInc=iinc;
% update newmark scheme parameters for the current increment
obj.updatenewmark(1);% important to use updatenewmark(theta=1) as of GN11 for static problem 03272019
% If Dt==DtO, then newmark parameters are the same, LHS and RHS are not to
% be changed ERROR: OBJ.DT SHOULD NOT BE TIMEINC(N+1)-TIMEINC(N)--060418
% FIXed on 060618
% As the enriched element Jacobian matrix is function to the crack
% displacement, crtLHS_UP has to be called every increment, changed on
% 11/29/18
obj.LinSysCrt.crtLHS_UP(obj.Newmark,iinc);     % update LHS for the current increment
% If newmark stays the same, the LHS and RHS stay the same for the linear system
% Cancel this filter on 07/10/2019 because the linear system may change
% size when the propagation is enabled.
obj.LinSysCrt.crtRHS_UP(obj.Newmark);          % update current external load vector;
obj.ExtLoadVec=obj.LinSysCrt.RHS;
% Zero U and P in linSysCrt to prepare for the new equation system in terms
% of the total unknown U(n+1) and P(n+1), Moved after crtLHS_UP on 11/30/18
% because crack opening update(which need U and P) is done within crtLHS_UP.
obj.switching(3);
% The next trial is very important to link step n+1 to step n
obj.LinSysCrt.upconf_trial(obj.Newmark);
obj.intforcer;                                     % Trial intforcer based on the first trial of u_sup(n+1) and p_sup(n+1)
% This residual load vector is the difference between external load vector and the trial internal load vector for the new time increment
obj.LinSysCrt.updateRHS(obj.ResLoadVec);
% Loop over equillibrium equations for a load increment
Ncut=0;                        % Number of times that increment cut procedure has been activated.
while obj.CutFlag==1
    while obj.ConvFlag==0
        % Check if load increment should be cut down
        % one way to switch the maxitr number between odd and even number.
        % 02/02/2021.
        if obj.IItr>=obj.MaxItr+Ncut || obj.DivFlag==1
            Ncut=Ncut+1;
            if Ncut>10
%                 cutratio=0.45;
%             elseif Ncut>10
                warning('The increment size may be too big, consider to use a small cut ratio.')
%                 cutratio=0.4;              
            end
            % cut the size of the current increment
            cutratio=0.4+0.2*rand(1);
            obj.autoincrem(1,cutratio,1.5,inclist);
            obj.updatenewmark(1);
            obj.LinSysCrt.crtLHS_UP(obj.Newmark,iinc);     % update LHS as newmark changed
            obj.LinSysCrt.crtRHS_UP(obj.Newmark);          % update external load vector;
            obj.ExtLoadVec=obj.LinSysCrt.RHS;
            obj.switching(3);
            obj.LinSysCrt.upconf_trial(obj.Newmark);
            obj.intforcer;
            obj.LinSysCrt.updateRHS(obj.ResLoadVec);
        end
        obj.IItr=obj.IItr+1;
        % Don't check stagechangeflag for the first iteration as it's not
        % accurate load update yet.
        if obj.IItr==1
            % Apply specified non-zero Dirichlet boundary condtion
            obj.LinSysCrt.appbound_v4; % applied on LHSnew
            % Use LinSysCrt.LHSnew (linear, based on total displacement)
            obj.LinSysCrt.solve(obj.Newmark,2); % solve the linearized equation system and update configuration.
            % Compute internal load force correponding to the external load vector at the current time step
            obj.intforcer;
        else
            % zero the Dirichlet boundary condition (04012019)
            obj.LinSysCrt.appbound_v5;  % applied on LHS
            % Use LinSysCrt.LHS (nonlinear, based on incremental displacement)
            obj.LinSysCrt.solve(obj.Newmark,1); % solve the linearized equation system and update configuration.
            % Compute internal load force correponding to the external load vector at the current time step
            stagechangeflag=obj.intforcer('stagecheck',1);
            % % when check stagechangeflag is true
            % Also return the stagechangeflag from cohesive traction update to
            % ensure that the updated cohesion traction does not exceed the
            % maximum cohesion 04152019.
            if stagechangeflag
                obj.ConvFlag=0;
                obj.DivFlag=1;
%                 fprintf('stagechange for cohesion updation, continue to inc cut\n')
                continue;
                % goto cutting time increments. Alternatively, one can cut the
                % load increments. It is just easier to cut time increments.
            end
        end
        obj.converger(psddofs); % NEED TO BE UPGRADED TO CHECK CONVERGENCE FOR U AND P SEPARATELY-- 05112018
        %         obj.ConvFlag=1;
    end
    % calculate the max and min of the pressure incrementand check if the
    % increment size should be increased or decreased
    obj.checksize_porepressure(stdpdofs,Ncut,psddofs,inclist);
    % Do crack growth checking only when the Cut is not needed
    if obj.CutFlag~=1
        % Check if the time increment is too large to lead unstable crack
        % growth prediction. If cutflag is true, then cut the current time
        % increment.  07112019
        if ~isempty(obj.LinSysCrt.EnrichItems)
            unstablegrowflag=false(1,length(obj.LinSysCrt.EnrichItems));
            for ienr=1:length(obj.LinSysCrt.EnrichItems)
                [unstablegrowflag(ienr),cutflag]=obj.LinSysCrt.EnrichItems(ienr).check_grow;
                if cutflag
                    obj.ConvFlag=0;
                    obj.DivFlag=1;
                    obj.CutFlag=1;
                    % BUG Fixing: when there are multiple enrichitems, we
                    % should break the loop but not continue to the next
                    % iteration. This may induce slow increment as only one
                    % cutflag has to delay other crack grow but it is
                    % required to have convergent solution.11/06/2020.
                    break; 
                end
            end
        end
    end
end
% % Update EnrichItems after the converged solution (May not be accurate
% enough). % IF NEEDED TO UPDATE ENRICHITEMS WITHIN THE ITERATION IS NOT
% CONFIRMED. (03102019) Attempt to change it to update enrichitems within
% iterations, suspended on 03132019 because much need to changed to
% accomodate it. obj.update_enrich; % Comprehensive method: postprocess
% cracks and grow cracks update_enrich is now a method of DomainPack.Domain
% as the update comes after the converged solution. There is no need to do
% it within iterating.m.11/27/20.

% adjust the increment size if unstablegrow is requested.
if any(unstablegrowflag)
    % minimal time step =inc*obj.Tottime, must be smaller than the step
    % size defined in the main function.
    minimalinc=1e-5;    
    allowedsteps=2;
    obj.autoincrem(3,0.25,1.5,inclist,minimalinc,allowedsteps);
end
% store current values of some parameters to the last converged value
obj.switching(1);
fprintf('No.%d increment finished within %d iterations\n',iinc,obj.IItr);
fprintf('The progress percentage is %f%% \n', obj.Increments(iinc)*100);
end





