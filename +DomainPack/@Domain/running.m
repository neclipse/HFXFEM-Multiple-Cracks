%% Method of Domain class to direct the whole solving process
function obj=running(obj,postdict,savemode,varargin)
%RUNNING Method of 'Domain' class
% Loop over equilibrium iterations for one load increment
% Assemble, solve Linear equation system and update
p=inputParser();
defaultsaveinterval=1;
defaultinclist=[];
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
% 10/30/20 There may be more than one encracks in the future, so that
% the following initialization would be carried out in a domain class.
% Another reason is that the initial enrich and update_enrich will be
% carried out in two steps: first we set_enrich and second we use
% subdomain to the NewElems for EnrichGauss and enrich these elems.
obj.initiate_enrich;
obj=obj.updatedofarray_enriched; % Bug-Issue #15. 12/04/2020
obj.LinSysCrt.initialRHS;
% Add try and catch statements to handle unexpected errors without saving
% the results to obj.Postprocess.07132019

while iinc<=obj.NewtonRaphson.NoInc                 % For loop is not proper for this kind loop whose loop number will change during the loop
    %         fprintf('Running at No. %d increment out of %d increments .\n',iinc,obj.NewtonRaphson.NoInc);
    % run into each inc with NewtonRaphson.iterating.
    % include both standard and enriched dofs described by Dirichlet boundary condition.
    % BUG-ISSUE 13: allpsddofs did not correctly include the
    % psdenrdofs. Reason is that obj.initiate_enrich did not return
    % obj. Remember Domain is a value class.
    allpsddofs=[obj.PsdDofs';obj.PsdEnrDofs];
    obj.NewtonRaphson.iterating(iinc,stdpdofs,allpsddofs,inclist);
    % Do postprocessing after convergence, then update the enrichitems.
    obj.update_enrich;
    %%- Store the converged value for postprocessing
    inc=obj.NewtonRaphson.Timeinc(iinc);            % the value of current time
    if savemode==1
        % iterating over one single increment
        obj=obj.storage(postdict,iinc,inc);
        ind=iinc;
    elseif savemode==2
        % store the results every specified inc
        if iinc==1
            obj=obj.storage(postdict,iinc,inc);
        elseif ~mod(iinc,saveinc)||iinc==obj.NewtonRaphson.NoInc
            ind=ind+1;
            obj=obj.storage(postdict,iinc,inc,ind);
        end
    elseif savemode==3
        % store the results at the specified increment
        if ismembertol(inc,timelist)
            obj=obj.storage(postdict,iinc,inc,ind);
            ind=ind+1;
        end
    end
    % Early termination due to cut through
    % This module should be updated because the false Isactive may not
    % be equivalent to cut through when there are multiple cracks.
    % 12/07/20.
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

if savemode==3
    postdict=postdict(1:ind-1);
else
    postdict=postdict(1:ind);
end
obj.Postprocess=postdict;
end
