function checksize_porepressure ( obj,stdpdofs,varargin )
%CHECKSIZE Control the size of next increment
% Author: Chang Huang, LSU. email: huangchang73@gmail.com
%    conisderations:
% 1. The variation of pore pressure
%% Input parser
p=inputParser;
addRequired(p,'Ncut');
addOptional(p,'psddofs',[]);
addOptional(p,'timelist',[]);
% addOptional(p,'unstablegrow',false);
parse(p,varargin{:});
Ncut=p.Results.Ncut;
psddofs=p.Results.psddofs;
timelist=p.Results.timelist;

%% Preparation based on pore pressure
P=obj.LinSysCrt.P;
PO=obj.LinSysCrt.PO;
% retract the pindex of the potential predofs for pressure
[~,ppsddofs]=mapping2(psddofs);
P(ppsddofs)=0;
PO(ppsddofs)=0;
Pt1=obj.LinSysCrt.Pt1;
Pt1(ppsddofs)=0;
% exclude the enriched pore pressure dofs, added on 08202019
P=P(1:stdpdofs);
PO=PO(1:stdpdofs);
Pt1=Pt1(1:stdpdofs);
% Find the maximum change in pore pressure in this step
maxpinc=max(P-PO);
minpinc=min(P-PO);
mpinc=[maxpinc,minpinc];
maxprate=max(abs(Pt1));
[absmaxpinc,~]=max(abs(mpinc)); 
%% try to enlarge the next increment
if obj.IInc>1 && obj.IInc< obj.NoInc
    if ~isempty(obj.Pincallowed)
        % Decrease the size of the next increment
        if absmaxpinc>obj.Pincallowed
            obj.DivFlag=1;
            obj.ConvFlag=0;
            % Goto Iterating: if obj.DivFlag==1, then cut increment
        else
            obj.CutFlag=0;
            % store max and min of the pressure increment
            obj.Pincmax=[obj.Pincmax;mpinc];
            obj.Pratemax=[obj.Pratemax;maxprate];
        end
    else
        obj.CutFlag=0;
        % store max and min of the pressure increment
        obj.Pincmax=[obj.Pincmax;mpinc];
        obj.Pratemax=[obj.Pratemax;maxprate];
    end
    if Ncut==0 && ~isempty(obj.Pinclimit)
        [absmaxpinc,~]=max(abs(mpinc));
        % Elongate the size of the next increment
        if absmaxpinc<obj.Pinclimit
            % Add a Test: if the last increment has crack growth 082019
            if ~isempty(obj.LinSysCrt.EnrichItems)
%                 newelems=[obj.LinSysCrt.EnrichItems{:}.NewElems];
%                 if isempty(newelems)
%                     obj.autoincrem(2,0.2,1.4,timelist);
%                 end
            else
                obj.autoincrem(2,0.2,1.5,timelist);
            end
            % Then break from the while obj.CutFlag==1 loop in iterating.
        end
    end
else
    obj.CutFlag=0;
    % store max and min of the pressure increment
    obj.Pincmax=[obj.Pincmax;mpinc];
    obj.Pratemax=[obj.Pratemax;maxprate];
end

end

