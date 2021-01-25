function [ gp,gw ] = linegauss( obj,id,cohesive,perforated,varargin )
% LINEGAUSS is a method of Enriched Elem_2d_UP, called by EnrCrackBody
% Use the seeds to create local gaussian points for line integral
% The number of line integral gaussian points are better to be at least 2 as 
% int(N'*N,dtau) is encountered. Also, the local line gaussian points ksi1 and
% ksi2 should be transformed to local elemental coordinates by using the line 
% shape function and the local element coordinates of the two intersections.
% It is also not wise to use too many line gaussian points to avoid
% spurious oscillations of traction  in the adherent zone (Ferte et al 2016).
% verified on 03292019
inp=inputParser;
addOptional(inp,'Alpha',pi/2);
addOptional(inp,'p',2);
parse(inp,varargin{:});
Alpha=inp.Results.Alpha;
p=inp.Results.p;
% global lcr dcr tmax lini;
% agl=assembleglobalinputs();
gbinp=obj.GaussPntDictM(1).GBINP;
lcr=gbinp.lcr;
dcr=gbinp.dcr;
tkrg=gbinp.tkrg;
lini=gbinp.lini;
tini=gbinp.threshold;
minaperture=gbinp.minaperture;
perfaperture=gbinp.perfaperture;
% use logical array id to relieve the single index id. 10/02/20
id=obj.Enrich==id;
localint=obj.LocalInt{id}; % elemental local coordinates of the intersections [xi,eta]
[ksi,lw]=GaussQuad(p);% parent segment [-1,1]
Nline=[1-ksi,1+ksi]/2;
gp=Nline*localint;
gw=lw;
nnodes=length(obj.NodList);
% calculate the normal and the tangent vector of the crack
% as of the straight line assumption of the crack segment within an
% element, these vectors are independent on the gaussian points
Nline_xi=[-1/2,1/2];
globalpnts=obj.GlobalInt{id};
nx=Nline_xi*(globalpnts(:,2));                  % y_xi
ny=-Nline_xi*(globalpnts(:,1));                 % -x_xi
mx=ny;                                          % -x_xi
my=-nx;                                         % -y_xi
ds=sqrt(nx^2+ny^2);
% THE DIRECTION OF NTAUD MAY MATTER 11/05/2018
% THE DIRECTION WAS MADE CONSISTENT IN "INTERSECTION.M" ON 11/06/18
% The direction of ntaud is consistent with the level set calculation in
% opengeo. confirmed on 12/16/20. 
%No, actually it is not consistent as the order of globalpnts changed when
%the order of opengeo.segments was from right to left. Fixed in
%opengeo.discretize, as for issue #32.
ntaud=[nx;ny];
mtaud=[mx;my];
Ntaud=ntaud/ds;   % unit normal vector
Mtaud=mtaud/ds;   % unit tangent vector

% NOTE the order of the gauss points in these enriched elements
%are not consistent with those in traditional elements (1,3,4,2)

gaussdictm(1,p)=FEPack.GaussPnt_Cohesive(); % generate an void object array
for igauss=1:p
    gaussdictm(igauss)=FEPack.GaussPnt_Cohesive(gp(igauss,1),gp(igauss,2),gw(igauss),nnodes,gbinp,perforated,Alpha);
    gaussdictm(igauss)=gaussdictm(igauss).preparing(obj.X,obj.Y,obj.EnrichNum);
    % Explicit linear softening
    if strcmp(cohesive,'unified')
        gaussdictm(igauss).TractionLaw=EnrichPack.TractionLaw_unified(lcr,dcr,tkrg,tini,lini,gbinp);
    elseif  strcmp(cohesive,'linear')
        gaussdictm(igauss).TractionLaw=EnrichPack.TractionLaw_linearsoftening(dcr,tkrg,lini,gbinp);
    elseif strcmp(cohesive,'bilinear')
        % Include hardening also 07052019
        gaussdictm(igauss).TractionLaw=EnrichPack.TractionLaw_bilinearsoftening(lcr,dcr,tkrg,lini,gbinp);
    end
    gaussdictm(igauss).Ntaud=Ntaud;
    gaussdictm(igauss).Mtaud=Mtaud;
    gaussdictm(igauss).Ds=ds;
    gaussdictm(igauss)=gaussdictm(igauss).initiate(minaperture,perfaperture);
%     gaussdictm(igauss)=gaussdictm(igauss).initialize; % 04042019 not finished
end
obj.LineGaussDict{id}=gaussdictm;
end

