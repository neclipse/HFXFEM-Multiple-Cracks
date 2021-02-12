function [unstablegrow,cutflag]=check_grow(obj,varargin)
% concrete method of EnrCrackBody
% update geometry and associated elements and nodes, selectively enrich new elements and nodes
% check the propagation rule using the component obj.Propagate (not realised)
% assign basic enrichment info to the selected elements and
% nodes
%% Lead obj.Mytips to lookahead
growflags=false(1,length(obj.Mytips));
unstablegrowflags=false(1,length(obj.Mytips));
cutflags=false(1,length(obj.Mytips));
unstablegrow=false;
cutflag=false;
for itip=1:length(obj.Mytips)
    if any(obj.Stdnodes(itip,:)) % the call of Stdnodes will check the smeared flag of tip element.
        % calculate the stress at tip
        obj.Mytips(itip).calstress_nonlocal;
        % look ahead to see if the crack shall propagate
        [growflags(itip),unstablegrowflags(itip),cutflags(itip)]=obj.Mytips(itip).lookahead;
    end
end
if any(cutflags)
%     disp('cut back the current time increment');
    cutflag=true;
    return;
end
if any(unstablegrowflags)
%     disp('cut the following increments for unstable growth');
    unstablegrow=true;      % allow grow but need reduce inc for the following increments
end
end

