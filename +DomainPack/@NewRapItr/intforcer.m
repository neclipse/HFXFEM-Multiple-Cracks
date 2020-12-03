function stagechangeflag=intforcer( obj, varargin)
%INTFORCER Method of Returner class
%   Loop over all elements to calculate and assemble the global interal
%   load vector, this method will call elemental method--ifstd
p=inputParser();
defaultstagecheck=0;
addRequired(p,'obj');
addOptional(p,'stagecheck',defaultstagecheck);
parse(p,obj,varargin{:});
obj=p.Results.obj;
stagecheck=p.Results.stagecheck;
edict=obj.LinSysCrt.ElemDict;
% IMPORTANT ERROR: LOADVEC SHOULD BE RESET TO ZERO FOR EVERY ITERATION SINCE 
% THE ELEMENT INTLOADVEC IS ALREADY ACCUMULATED LOAD, INSTEAD OF LOAD CHANGE
loadvec=zeros(obj.Dim,1);
calstress=false;
for ielem=1:length(edict)
    elem=edict(ielem);
    if any(elem.Enrich)     % Enriched elements
%         ienrich= find(elem.Enrich,1);
        % At the preliminary stage, there will be at most one enrich item in
        % one element, so this works. Later on, when there can be more enrich
        % items in one element, it has to be changed. 11/22/2018
        % To change 091120
        % chaged on 10/22/20, now the load update is elemental. No need to
        % loop over every enrichitem.
        if stagecheck
            [load,stagechangeflag]=elem.ifstd_enriched(obj.Newmark,stagecheck,calstress);
            if stagechangeflag
                return;  % return without proceedeing to the following
            end
        else
            load=elem.ifstd_enriched(obj.Newmark,stagecheck,calstress);
        end
        array=elem.JacobianMat.LocarrayAll;
    else                    % Standard elements
        stagechangeflag=false;
        load=elem.ifstd(obj.Newmark,calstress);
        array=edict(ielem).Locarray;
    end    
    loadvec(array)=loadvec(array)+load;
end
obj.IntLoadVec=loadvec;
end
