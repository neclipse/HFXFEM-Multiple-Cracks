function  postprocess( obj, dt )
%POSTPROCESS of the EnrCrackBody is to obtain the practical fracture
%information from the converged solution after each increment.
% calculate tip stress
%   Loop through the INTELEM and the linegauss
%% Clean slate
agl=obj.Elemdict(1).GaussPntDictM(1).GBINP;
inipore=agl.inipore;
% skin=agl.skin;
numle=length(obj.INTELEM);
numlg=length(obj.INTELEM(1).LineGaussDict{obj.Id});
numl=numle*numlg;
% global coordinates (x,y) of the line gaussian integration points
intpoints=zeros(numl,2);
uplus=zeros(numl,2);
uminus=uplus;
aperture=zeros(numl,1);
pfrack=zeros(numl,1);
% %cohesion traction (MPa) is not Jacob.F_coh, as F_coh is the cohesion Force (N)
% % integrated along the discontinuity, traction*area=Force
ctraction=zeros(numl,2);    % cohesion traction  [tx1,ty1; tx2, ty2;... ];
crackvf=zeros(numl,1);      % integration fraction for every
% leakvfstd=zeros(numle,obj.INTELEM(1).NoNodes);       % distribution of leak-off FLUX from the standdofs to each node
% leakvfenr=zeros(numle,obj.INTELEM(1).NoNodes);       % distribution of leak-off FLUX from the Enricheddofs to each node
iint=1;
for in=1:length(obj.ENNODE)
    obj.ENNODE(in).Leakoff= 0;
end
for iE=1:length(obj.INTELEM)
    %% Withdraw data
    elem=obj.INTELEM(iE);
    Jacob=elem.JacobianMat;
    % ua is the current displacement discontinuity of all nodes for the
    % specified (id) enrichement item.
    % Standard Dofs
    us=elem.Un2i;
    ps=elem.Pn2i;      % elemental iterative total pore pressure vector, p_sup(n+1)_sub(i+1)
    pds=elem.Pn2t1;
    % due=Jacob.dUn2iEnr;
    ue=Jacob.Un2iEnr;
    pe=Jacob.Pn2iEnr;
    ve=Jacob.Un2t1Enr;
    pde=Jacob.Pn2t1Enr;
    % Enriched Dofs, zero initial values if they are newly enriched, added
    % on 07132019
    % 11/06/20 The NewElems are made empty by NewRapItr.update_enrich
    % because the postprocessing comes before real crack growth.
%     if ~isempty(obj.NewElems)
%         if any(obj.NewElems==elem.Ind)
%             %  due=zeros(size(us));
%             ue=zeros(size(us));
%             pe=zeros(size(ps));
%             ve=zeros(size(us));
%             pde=zeros(size(pds));
%             
%         end
%     end
    %% Leakoff calculation
    % Integrate the leak off volume, from two parts: standard and enriched
    % NOTE HERE JACOB MATRICES ARE UPDATED FROM THE UPDATED CRACK
    % OPENING. (03132019)
    % It is confirmed that qint has such a sign convention: positive
    % for flow into the domain from the crack; negative for flow
    % from the domain to the crack.(04232019)
%     if ~isempty(obj.NewElems)
%         if any(obj.NewElems==elem.Ind)
%             qintps=zeros(size(ps));
%             qintpe=zeros(size(ps));
%         else
%             qintps=Jacob.Hintpsps*ps+Jacob.Hintpspe*pe+Jacob.Sintpsps*pds+Jacob.Sintpspe*pde-Jacob.Qintueps'*ve;
%             qintpe=Jacob.Hintpeps*ps+Jacob.Hintpepe*pe+Jacob.Sintpeps*pds+Jacob.Sintpepe*pde-Jacob.Qintuepe'*ve;
%         end
%     else
        qintps=Jacob.Hintpsps*ps+Jacob.Hintpspe*pe+Jacob.Sintpsps*pds+Jacob.Sintpspe*pde-Jacob.Qintueps'*ve;
        qintpe=Jacob.Hintpeps*ps+Jacob.Hintpepe*pe+Jacob.Sintpeps*pds+Jacob.Sintpepe*pde-Jacob.Qintuepe'*ve;
%     end
    %    leakvfstd(iE,:)=transpose(qintps);      % The leak-off FLUX distributed to the standard dofs of the nodes within iE element
    %    leakvfenr(iE,:)=transpose(qintpe);      % The leak-off FLUX distributed to the enriched dofs of the nodes within iE element
    % 04292019, need subtract the node.Leakoff calculated from the
    % cal_qextenr.
    % Bug: should sum on the Leakoff but not CInjection, fixed on 05092019
    % Bug: should zero NodDict.Leakoff in every step. 05302019, changed back
    % to old method for efficiency.
    for in=1:elem.NoNodes
        elem.NodDict(in).Leakoff= elem.NodDict(in).Leakoff+qintps(in)+qintpe(in);
        elem.NodDict(in).ALeakoff= elem.NodDict(in).ALeakoff+(qintps(in)+qintpe(in))*dt;
    end
    
    %% Calculation at the linegaussian points along the crack curve
    linegauss=elem.LineGaussDict{elem.Enrich==obj.Id};
    for ig=1:length(linegauss)
        lg=linegauss(ig);
        intpoints(iint,:)=[lg.X,lg.Y];
        % update the crack opening using enriched udofs at the begining (03222019)
        lg=lg.updatecrackopening(us,ue);
        % 08232019: Add +obj.IniCrackDisp to Uplus so that the final
        % crackdisp and the crackopening also accounts for the initial
        % setting crack aperture only in postprocess as gausscohesive.uplus
        % may be used elsewhere.
        uplus(iint,:)=lg.Uplus'+1/2*lg.IniCrackDisp';
        uminus(iint,:)=lg.Uminus'-1/2*lg.IniCrackDisp';
        aperture(iint)=lg.CrackOpening;
        pf=lg.Np*ps+lg.Npenr*pe+inipore;      % should use the enriched interpolation
        pfrack(iint)=pf;
        %        % use the incremental form to update traction vector;
        %        lg.Traction=lg.Traction+lg.Tangent_coh*lg.Nu*due;
        ctraction(iint,:)=transpose(lg.Traction);
        % Integrate the fracture volume (fraction on the current gauss point)
        crackvf(iint)=lg.H*(lg.CrackOpening-lg.MinCrackOpening)*lg.Ds;
        % update the crack opening using enriched udofs at the end (03142019)
        % NEED RETURN the lg to elem.LineGaussDict{obj.Id} 03222019
        elem.LineGaussDict{obj.Id}(ig)=lg;
        iint=iint+1;
    end
    
end
%% Sort and Store the results
% cinjection is calculated from cal_qextenr.
cinjection=[obj.ENNODE.CInjection];       % positive is injection, the sign is already flipped in cal_qextenr
allleakvf=[obj.ENNODE.Leakoff];
effectivelkf=sum(allleakvf)+sum(cinjection);
% leakvf=allleakvf(cinjection==0);
% leakvf=leakvfstd+leakvfenr;             % The combined leakoff volume spreaded to the discrete nodes of the element, (80% percent sure 03182019)
% LEAKOFF FLUX INTEGRATED WITH TIME GIVES THE LEAKOFF VOLUME (NOT FINISHED,04022019)
% reorder the integration points and other practical information
tip=obj.Mygeo.Tips(2,:);
temp=intpoints-repmat(tip,numl,1);
dist=sqrt(temp(:,1).^2+temp(:,2).^2);
[~,I]=sort(dist);
obj.IntPoints=intpoints(I,:);
obj.Uplus=uplus(I,:);
obj.Uminus=uminus(I,:);
obj.Aperture=aperture(I);
obj.Pfrack=pfrack(I);
obj.CTraction=ctraction(I,:);
obj.CrackVf=crackvf(I);
obj.CrackVolume=sum(crackvf);
% monitoring leak-off FLUX, to recover Carter's leakoff coefficient
% mesh01082020
% obj.LeakoffFlux=[obj.Nodedict([574,584,594,604]).Leakoff];% m^2/s
% mesh 01122020
% obj.LeakoffFlux=[obj.Nodedict([9333,9317,9297,9277,9257,9237,9217,9207]).Leakoff];% m^2/s
% mesh 040320 Dontsov_medium
% obj.LeakoffFlux=[obj.Nodedict([445,4836,1615,1635,1645,1705,1755,1865,1915,1945]).Leakoff];% m^2/s
obj.LeakoffVolume=obj.LeakoffVolume+effectivelkf*dt; % approximated total leak-off volume
obj.Length=obj.Mygeo.Length;
obj.CMOD=obj.Aperture(1);
obj.CMP=obj.Pfrack(1);
end

