function peakTangent=givepeaktangent(obj)
scr=obj.CriticalDisp;
dcr=obj.PeakTraction;
%Non-dimension shear and normal components
Tss=-scr/dcr;        % an arbitray large positive number, changed to zero
Tsn=0;
Tnn=-scr/dcr;   % initialize for newly generated cohesive segments
peakTangent=[Tss,Tsn;Tsn,Tnn];
end