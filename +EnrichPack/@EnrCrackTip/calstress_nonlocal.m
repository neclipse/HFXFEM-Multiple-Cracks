function calstress_nonlocal(obj)
% calculate the stress at the tip using nonlocal gaussian weights
% (Wells and Sluys, 2001)
% 1. Roughly locate the elements within the interaction radius
if obj.Isactive
    dis=sqrt((obj.Mesh.VX-obj.Xct).^2+(obj.Mesh.VY-obj.Yct).^2);
    foundnodes=find(dis<obj.Radius);
    possibleelems=zeros(200,1);
    iloc1=1;
    for in=1:length(foundnodes)
        foundelems=obj.Mesh.nodetoelem(foundnodes(in));
        iloc2=iloc1+length(foundelems)-1;
        possibleelems(iloc1:iloc2)=foundelems;
        iloc1=iloc2+1;
    end
    possibleelems=possibleelems(1:iloc2);
    % Exclude the repetition
    elempool=unique(possibleelems);
    % 2. locate the gaussian points within the interaction radius
    % (ahead of the crack tip)
    % BUG: THE GAUSSPNDICTM DOES NOT HAVE STRESS VALUES CALCULATED AT
    % ENRICHED ELEMENTS.
    % Quick Fix:delete the enriched elements in advance 08112019
    %     elempool=setdiff(elempool_temp,obj.Mygeo.Intelements);
    
    % Real fix: call ifstd within the method for all elems in elempool_temp
    % 09222019, replaced by calstress, 12/22/2020
    storage=false;
    gausspool=[];
    for ie=1:length(elempool)
        if any(obj.Elemdict(elempool(ie)).Enrich)     % Enriched elements
            % HERE ONLY ONE ENRCRACK IS ASSUMED. TO CHANGE. 09/18/20
            % Fixed and ienrich is obsolete. 11/20
%             ienrich= find(obj.Elemdict(elempool(ie)).Enrich,1);
            obj.Elemdict(elempool(ie)).calstress_enriched(storage);
            gausspool=[gausspool,obj.Elemdict(elempool(ie)).EnrichGauss];
        else
            obj.Elemdict(elempool(ie)).calstress(storage);
            gausspool=[gausspool,obj.Elemdict(elempool(ie)).GaussPntDictM];
        end
    end
%     gausspool=[obj.Elemdict(elempool).GaussPntDictM];
    X=[gausspool.X];
    Y=[gausspool.Y];
    [xyl,r]=obj.callocal(X,Y);
    % a. test if the gaussian point is within the radius
    inside=find(r<obj.Radius);
    % b. test if the gaussian point is ahead of the crack tip
    ahead=find(xyl(1,:)>0);       % one column one point
    goodpointsind=intersect(inside,ahead);
%     goodpointsind=inside;
    goodpoints=gausspool(goodpointsind);
    goodr=r(goodpointsind);
%     goodpoints=gausspool;
%     goodr=r;
    % 3. retrieve the effective stress data at these points
    stressp=[goodpoints.Stressp]; % one column one point
    % 4. get the weighted average
    l=obj.Radius;
    w=1/((2*pi)^1.5*l^3)*exp(-goodr.^2/(2*l^2));
    weights=sum(w);
    sx=stressp(1,:)*w'/weights;
    sy=stressp(2,:)*w'/weights;
    sxy=stressp(3,:)*w'/weights;
    sz=stressp(4,:)*w'/weights;
    obj.Stressp=[sx;sy;sxy;sz];
end
end

