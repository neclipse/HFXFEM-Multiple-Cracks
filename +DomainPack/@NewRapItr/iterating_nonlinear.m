function  iterating_nonlinear( obj, iinc, varargin)
%Iterating Method of NewtonRaphson class
% Loop over equilibrium iterations for one increment, modified
% Newton-Raphson method as the tangent matrix is updated per increment.
% Assemble, solve Linear equation system and update
if nargin>3
    psddofs=varargin{1};
    timelist=varargin{2};
elseif nargin>2
    psddofs=varargin{1};
    timelist=[];
else
    psddofs=[];
    timelist=[];
end
obj.IInc=iinc;

% update newmark scheme parameters for the current increment
obj.updatenewmark(1);% important to use updatenewmark(1) as of GN11 for static problem 03272019
% If Dt==DtO, then newmark parameters are the same, LHS and RHS are not to
% be changed ERROR: OBJ.DT SHOULD NOT BE TIMEINC(N+1)-TIMEINC(N)--060418
% FIXed on 060618 
% As the enriched element Jacobian matrix is function to the crack
% displacement, crtLHS_UP has to be called every increment, changed on
% 11/29/18
obj.LinSysCrt.crtLHS_UP(obj.Newmark,iinc);     % update LHS for the current increment
 % If newmark stays the same, the LHS and RHS stay the same for the linear system 
if abs(obj.Dt-obj.DtO)>1e-6                      
obj.LinSysCrt.crtRHS_UP(obj.Newmark);          % update current external load vector;
obj.ExtLoadVec=obj.LinSysCrt.RHS;
end
% Zero U and P in linSysCrt to prepare for the new equation system in terms
% of the total unknown U(n+1) and P(n+1), Moved after crtLHS_UP on 11/30/18
% because crack opening update(which need U and P) is done within crtLHS_UP.
obj.switching(3);
% The next trial is very important to link step n+1 to step n
obj.LinSysCrt.upconf_trial(obj.Newmark);           
obj.intforcer;                                     % Trial intforcer based on the first trial of u_sup(n+1) and p_sup(n+1)
% This residual load vector is the difference between external load vector and the trial internal load vector for the new time increment
obj.LinSysCrt.RHS=obj.ResLoadVec;                      
% Loop over equillibrium equations for a load increment
Ncut=0;                        % Number of times that increment cut procedure has been activated.
while obj.CutFlag==1
    while obj.ConvFlag==0
        % Check if load increment should be cut down
        if obj.IItr>=obj.MaxItr || obj.DivFlag==1
            Ncut=Ncut+1;
            if Ncut>5
                error('The simulation aborted because the increment size is not set to get a convergence.')
            end
            % cut the size of the current increment
            obj.autoincrem(1,0.2,1.5,timelist);
            obj.updatenewmark(0.99);
            obj.LinSysCrt.crtLHS_UP(obj.Newmark,iinc);     % update LHS
            obj.LinSysCrt.crtRHS_UP(obj.Newmark);          % update external load vector;
            obj.ExtLoadVec=obj.LinSysCrt.RHS;
            obj.switching(3);
            obj.LinSysCrt.upconf_trial(obj.Newmark);
            obj.intforcer;                                
            obj.LinSysCrt.RHS=obj.ResLoadVec;
        end
        obj.IItr=obj.IItr+1;
        if obj.IItr==1
            % Apply specified non-zero Dirichlet boundary condtion
            obj.LinSysCrt.appbound_v4; % applied on LHSnew
            % Use LinSysCrt.LHSnew (linear, based on total displacement)
            obj.LinSysCrt.solve(obj.Newmark,2); % solve the linearized equation system and update configuration.
        else
            % zero the Dirichlet boundary condition (04012019)
            obj.LinSysCrt.appbound_v5;  % applied on LHS
            % Use LinSysCrt.LHS (nonlinear, based on incremental displacement)
            obj.LinSysCrt.solve(obj.Newmark,1); % solve the linearized equation system and update configuration.
        end
        % Compute internal load force correponding to the external load vector at the current time step
        obj.intforcer; % Need to update T_coh after the solution 0322, (NOT FIXED)
        obj.converger(psddofs); % NEED TO BE UPGRADED TO CHECK CONVERGENCE FOR U AND P SEPARATELY-- 05112018
%         obj.ConvFlag=1;
    end
    % calculate the max and min of the pressure incrementand check if the
    % increment size should be increased or decreased
    obj.checksize(psddofs,timelist);
end
% % Update EnrichItems after the converged solution (May not be accurate enough). 
% % IF NEEDED TO UPDATE ENRICHITEMS WITHIN THE ITERATION IS NOT CONFIRMED. (03102019)
% Attempt to change it to update enrichitems within iterations, suspended
% on 03132019 because much need to changed to accomodate it.
if ~isempty(obj.LinSysCrt.EnrichItems)
    for ienr=1:length(obj.LinSysCrt.EnrichItems)
        obj.LinSysCrt.EnrichItems{ienr}.postprocess(obj.Dt);      % obtain crack aperture and other practical information.
    end
end
% store current values of some parameters to the last converged value
obj.switching(1);
fprintf('No.%d increment finished within %d iterations\n',iinc,obj.IItr);
end





