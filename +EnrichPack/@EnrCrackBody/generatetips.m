function generatetips( obj )
%GENERATETIPS Construct EnrCrackTip objects and assign to obj.Mytips;
% 1. Check Mygeo to detect the real tips
% 2. Construct the encracktip objects corresponding to the real tips;
Rtips=obj.Mygeo.Rtips;
% test if there is a real tip
if ~isempty(Rtips)
    mytips=[];
    % loop over the Rtips and construct EnrCrackTip object
    for i=1:length(Rtips)           
        itip=Rtips(i);              % index of the real tip
        mytip=EnrichPack.EnrCrackTip('cracktip',obj.Elemdict,obj.Nodedict,obj.Mygeo,itip);
        mytip.Id=obj.Id; % tip has the same Id as its crackbody
        mytips=[mytips,mytip];
    end
    obj.Mytips=mytips;
end

end

