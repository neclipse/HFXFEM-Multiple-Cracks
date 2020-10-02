function [IntLoadAll,stagechangeflag]=ifstd_enriched( obj,newmark,id,stagecheck,calstress)
% internal force vector calculation and the stress update
% global skin;
% skin=gbinp.skin;
	% Enrichement screening
% 	if obj.Enrich(id)==0
% 		error('The element does not interact with the specified enrichement item\n')
% 	end 
    % use logical array id to relieve the single index id. 10/02/20
    id=obj.Enrich==id;
    
	% Standard Dofs
    us=obj.Un2i;      % elemental iterative total displacement vector, u_sup(n+1)_sub(i+1)
    ps=obj.Pn2i;      % elemental iterative total pore pressure vector, p_sup(n+1)_sub(i+1)
    vs=obj.Un2t1;
    as=obj.Un2t2;
    pds=obj.Pn2t1;
	% Enriched Dofs
	Jacob=obj.JacobianMatDict(id);
    Due=Jacob.DUn2Enr;
	ue=Jacob.Un2iEnr;
	pe=Jacob.Pn2iEnr;   
	ve=Jacob.Un2t1Enr;
	ae=Jacob.Un2t2Enr;
	pde=Jacob.Pn2t1Enr;
	% Dynamic parameter: Newmark
    a1=newmark.a1;
    %% update stress at GaussPntDictM
    %  called by encracktip.calstress_nonlocal and elem.calstress in this
    %  version and the
    % newmark value is arbitrarily set as it is not used in matsu.
    if calstress
        for ig=1:length(obj.EnrichGaussDict{id})
            % update stresses
            obj.EnrichGaussDict{id}(ig)=obj.EnrichGaussDict{id}(ig).matsu_enriched(us,ue,ps,pe);
            % calculate the maximum principal stress and its orientation
            % [psmax,gtheta]=obj.EnrichGaussDict{id}.psu;
        end
        for ig=1:length(obj.LineGaussDict{id})
            % update stresses at the discontinuity
            obj.LineGaussDict{id}(ig)=obj.LineGaussDict{id}(ig).matsu_enriched(us,ue,ps,pe);
            % calculate the maximum principal stress and its orientation
            % [psmax,gtheta]=obj.EnrichGaussDict{id}.psu;
        end
        return;
        % return because the newmark value is arbitrarily set
    end
    %% update cohesive traction at linegauss_cohesive
    F_coh_i=zeros(size(Due));
    for ig=1:length(obj.LineGaussDict{id})
        lg= obj.LineGaussDict{id}(ig);   % VALUE CLASS, A HARD COPY
        if stagecheck
            % Do not return lg to lg
            [traction,stagechangeflag]=lg.matctu(us,ue,Due);
            if stagechangeflag
                IntLoadAll=[];
                return;
            end
        else
            %for the first iteration, the traction can inherit from the
            %last increment. 07272019. (Not 100% sure, changed back after trial)
            traction=lg.matctu(us,ue,Due);
%             traction=lg.TractionO;
            % For the new increment, do not assign traction to lg.Traction
            % yet. (does not matter now, because lg.Traction is updated based
            % on lg.TractionO and Due inside lg.matctu).
        end
        lg.Traction=traction;
        F_coh_i=F_coh_i+lg.H*lg.Nu'*lg.Traction*lg.Ds;
        % Need to add these for "value" class gauss
        obj.LineGaussDict{id}(ig)=lg;
    end
    %% Interal load vector for the element
    % calculate internal loadus, -a1*(Musus*as+Musue*ae+Kusus*us+Kusue*ue-Qusps*ps-Quspe*pe)
    loadus=-a1*(Jacob.Musus*as+Jacob.Musue*ae+Jacob.Kusus*us+Jacob.Kusue*ue-Jacob.Qusps*ps-Jacob.Quspe*pe);
    % Incremental update because Tangent_coh is nonlinear and updated once an increment
    %% 
    % calculate internal loadue, -a1*(Mueus*as+Mueue*ae+Kueus*us+Kueue*ue-Queps*ps-Quepe*pe+f_coh-f_int)
    % Incremental update as Jacob.Kc is nonlinear and updated once an increment
    % WRONG UPDATION FOR UNLOADING.(NOT FIXED) 03222019, Fixed on 04032019
    % by using incremental Due instead.
    % On 07/07/2019, the F_coh_old is not used at least for the first
    % increment when F_coh_old is actually empty.
%     if ~isempty(Jacob.F_coh_old)
%         Jacob.F_coh_i=Jacob.F_coh_old+Jacob.Kc*Due;
%     else
        Jacob.F_coh_i=F_coh_i;
%     end
    % To enforce the zero traction after complete separation
    if any(any(Jacob.Kc))==false % if all zero
        Jacob.F_coh_i=zeros(size(Due));
    end
	%f_int=Jacob.Qintueps*ps+Jacob.Qintuepe*pe;
	loadue=-a1*(Jacob.Mueus*as+Jacob.Mueue*ae+Jacob.Kueus*us...
	+Jacob.Kueue*ue-Jacob.Queps*ps-Jacob.Quepe*pe+Jacob.F_coh_i...
	-Jacob.Qintueps*ps-Jacob.Qintuepe*pe);
    % calculate internal loadp,
    % Qusps'*vs+Quspe'*ve+Hpsps*ps+Hpspe*pe+Spsps*pds+Spspe*pde-qintps; %
    % IMPORTANT ERROR 1: The minus sign is already included in the matrices
    % of Jacob, except Qintueps and Qintuepe. 04022019
	qintps=Jacob.Hintpsps*ps+Jacob.Hintpspe*pe+Jacob.Sintpsps*pds+Jacob.Sintpspe*pde-Jacob.Qintueps'*ve;
    % IMPORTANT ERROR 2: JACOB.QUSPE'*VE WAS A TYPO, CORRECTED TO
    % JACOB.QUEPS. THIS SHOULD ANSWER WHY LARGE ERROR ONLY OCCURS AT
    % STANDARD PDOFS WITHIN THE ENRICHED ELEMENT. (04022019)
    loadps=Jacob.Qusps'*vs+Jacob.Queps'*ve+Jacob.Hpsps*ps+Jacob.Hpspe*pe+Jacob.Spsps*pds+Jacob.Spspe*pde-qintps;
	% calculate internal loadpe, 
    qintpe=Jacob.Hintpeps*ps+Jacob.Hintpepe*pe+Jacob.Sintpeps*pds+Jacob.Sintpepe*pde-Jacob.Qintuepe'*ve;
    loadpe=Jacob.Quspe'*vs+Jacob.Quepe'*ve+Jacob.Hpeps*ps+Jacob.Hpepe*pe+Jacob.Speps*pds+Jacob.Spepe*pde-qintpe;
    %% need rearrange loadu and loadp to adapt to the data structure [fx1,fy1,fp1,fx2,fy2,fp2...,fx1enr,fy1enr,fp1enr...]
    % The approach is different from the rearrangement adopted in
    % crtLHS_UP_Uni, which is actually the opposite. The whole point is to
    % make sure the the local dof corresponds to its global
    % index.
    loadstd=zeros(length(Jacob.Locarray),1);
    loadstd(1:3:end-2)=loadus(1:2:end-1);       % extract from [fx1,fy1,fx2,fy2,...]
    loadstd(2:3:end-1)=loadus(2:2:end);
    loadstd(3:3:end)=loadps(:);                 % extract from [fp1,fp2,...]
    loadenr=zeros(length(Jacob.LocarrayEnr),1);
    loadenr(1:3:end-2)=loadue(1:2:end-1);
    loadenr(2:3:end-1)=loadue(2:2:end);
    loadenr(3:3:end)=loadpe(:);
    IntLoadAll=[loadstd;loadenr];            % comprehensive internal load vector corresponding to all dofs, [std;enr]
    Jacob.IntLoadVecAll=IntLoadAll;
    % [x1,y1,p1,x2, y2,p2...,x1enr,y1enr,p1enr...] is consistent with the
    % IntLoadAll
    Jacob.LocarrayAll=[Jacob.Locarray,Jacob.LocarrayEnr]; 
end

