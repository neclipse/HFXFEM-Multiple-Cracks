%% Method of Domain class to direct the whole solving process
function obj=running(obj,postdict,savemode,varargin)
%RUNNING Method of 'Domain' class
% Loop over equilibrium iterations for one load increment
% Assemble, solve Linear equation system and update
p=inputParser();
defaultsaveinterval=1;
defaultinclist=[0.01,0.1,1];
checkinterval = @(x) isinteger(x) && x<obj.NewtonRaphson.NoInc;
% addRequired(p,'savemode');
addOptional(p,'interval',defaultsaveinterval,checkinterval);
addOptional(p,'inclist',defaultinclist);
parse(p,varargin{:});
saveinc=p.Results.interval;
inclist=p.Results.inclist;
timelist=inclist*obj.NewtonRaphson.Tottime;
stdpdofs=obj.NoPDofs;
% maxinc=obj.NewtonRaphson.MaxInc;
% postdict(1,maxinc)=ToolPack.Postprocessor();
% obj.Postprocess=postdict;
% clear postdict;
%% Loop over Newton_Raphson iteration increments
ind=1;
iinc=1;
% obj.LinSysCrt.initialize; % initialize the dofs and calculate the initial external load vector due to inital stress;
obj.LinSysCrt.initialRHS;
% Add try and catch statements to handle unexpected errors without saving
% the results to obj.Postprocess.07132019
try
    while iinc<=obj.NewtonRaphson.NoInc                 % For loop is not proper for this kind loop whose loop number will change during the loop
%         fprintf('Running at No. %d increment out of %d increments .\n',iinc,obj.NewtonRaphson.NoInc);
        %%- Store the converged value for postprocessing
        inc=obj.NewtonRaphson.Timeinc(iinc);            % the value of current time
        allpsddofs=[obj.PsdDofs';obj.PsdEnrDofs];        % include both standard and enriched dofs described by Dirichlet boundary condition.
        if savemode==1
            % iterating over one single increment
            obj.NewtonRaphson.iterating(iinc,stdpdofs,allpsddofs);
            obj=obj.storage(postdict,iinc,inc);
            ind=iinc;
        elseif savemode==2
            % store the results every specified inc
            obj.NewtonRaphson.iterating(iinc,stdpdofs,allpsddofs);
            if iinc==1
                obj=obj.storage(postdict,iinc,inc);
            elseif ~mod(iinc,saveinc)||iinc==obj.NewtonRaphson.NoInc
                ind=ind+1;
                obj=obj.storage(postdict,iinc,inc,ind);
            end
        elseif savemode==3
            % store the results at the specified increment
            obj.NewtonRaphson.iterating(iinc,stdpdofs,allpsddofs,inclist);
            if ismembertol(inc,timelist)
                obj=obj.storage(postdict,iinc,inc,ind);
                ind=ind+1;
            end
        end
        if ~isempty(obj.EnrichItems)
            if ~all([obj.EnrichItems.Isactive])
                if obj.Postprocess(ind).IInc~=iinc
                    ind=ind+1;
                    obj=obj.storage(postdict,iinc,inc,ind);
                end
                fprintf('The simulation is ended early at %f seconds as crack cut through the domain\n', inc);
                break;
            end
            obj=obj.updatedofarray_enriched; % added on 06282019
            obj.LinSysCrt.initialRHS;
        end
        iinc=iinc+1;
    end
catch ME
    warning('Some unexpected error occured and the program exit before completion');
    % For debugging use, it is better to rethrow the error. Later on,
    % delete rethrow (ME).
    rethrow(ME);
end
if savemode==3
    postdict=postdict(1:ind-1);
else
    postdict=postdict(1:ind);
end
obj.Postprocess=postdict;
end
