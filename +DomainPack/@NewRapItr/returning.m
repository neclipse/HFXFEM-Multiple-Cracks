function  returning( obj )
%RETURNING Method of Retuner class
% Loop over equilibrium iterations for one load increment
% Assemble, solve Linear equation system and update 
% Evaluate external load vector

% Loop over Newton_Raphson iteration increments
obj.LinSysCrt.initialRHS(obj.Newmark);
iinc=1;
while iinc<obj.NoInc
    obj.IInc=iinc;
    % update newmark scheme parameters for the current increment
    obj.updatenewmark;
    % reset some parameters for converged problem
    obj.switching(2);
    obj.LinSysCrt.crtRHS_UP(obj.Newmark);          % update external load vector;
    obj.ExtLoadVec=obj.LinSysCrt.RHS;
    obj.LinSysCrt.upconf;                          % update the velocity and other solution dependent variable for the current increment
    obj.intforcer;
    obj.LinSysCrt.RHS=obj.ResLoadVec;
    obj.LinSysCrt.crtLHS(obj.Newmark,iinc);             % update LHS
% Loop over equillibrium equations for a load increment
    Ncut=0;                        % Number of times that increment cut procedure has been activated. 
    while obj.ConvFlag==0             
        % Check if load increment should be cut down
        if obj.IItr>=obj.MaxItr || obj.DivFlag==1
            Ncut=Ncut+1;
            if Ncut>5
                error('The simulation aborted because the increment size is not set to get a convergence.')
            end
            obj.IncCutFlag=1; 
            % cut the size of the current increment 
            obj.autoincrem(1);
            obj.updatenewmark;
            obj.switching(2);
            obj.LinSysCrt.crtRHS_UP(obj.Newmark);          % update external load vector;
            obj.ExtLoadVec=obj.LinSysCrt.RHS;
            obj.LinSysCrt.upconf;
            obj.intforcer;
            obj.LinSysCrt.RHS=obj.ResLoadVec;
            obj.LinSysCrt.crtLHS(obj.Newmark,iinc);             % update LHS
        end
        obj.IItr=obj.IItr+1;
        obj.LinSysCrt.solve(obj.Newmark);    
        % Compute global internal load force vector
        obj.intforcer;
        obj.converger;       
    end
    % calculate the max and min of the pressure increment
    maxpinc=max(obj.LinSysCrt.P-obj.LinSysCrt.PO);
    minpinc=min(obj.LinSysCrt.P-obj.LinSysCrt.PO);
    mpinc=[maxpinc,minpinc];
    maxprate=max(abs(obj.LinSysCrt.Pt1));
    [absmaxpinc,~]=max(abs(mpinc));
    % store max and min of the pressure increment
    obj.Pincmax=[obj.Pincmax;mpinc];        
    % store current values of some parameters to the last converged value
    obj.switching(1);
    % Elongate the size of the next increment
    if absmaxpinc<obj.Pincallowed && maxprate<obj.Prallowed && Ncut==0
       obj.autoincrem(2); 
    end
    iinc=iinc+1;
end

end



