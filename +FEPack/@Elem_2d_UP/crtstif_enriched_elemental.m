function crtstif_enriched_elemental( obj, newmark, gbinp, blending)
%CRTST method of elem_2d_EP to calculate the stiffness matrix using gauss
%point class
%   This function is slightly different from the calstifmatrix in class
%   elem_2d due to the new member GaussPntDict in the elem_2d_EP
%   Each gauss point object calculates the consistent tangent operator,D
%   matrix. This function is more like a assembler.
% global  Delastic Density Biot_alpha kmat Cstar Kf muf skin;
Delastic=gbinp.Delastic;
Density=gbinp.Density;
Biot_alpha=gbinp.Biot_alpha;
if blending
    kmat=gbinp.kmat_crack;
    Cstar=gbinp.Cstar_crack;
else
    % use smaller kmat for enriched elements to better control the fracture
    % flow within the crack plane, it may be reversed to domain kmat when
    % viscosity is relatively high.
    kmat=gbinp.kmat_crack;      
    Cstar=gbinp.Cstar;          % use consistent domain Biot_mod
end
muf=gbinp.muf;
k=obj.EnrichNum;                % number of current enrichitems in the element
nnodes=obj.NoNodes;
% 10/09/20: This scheme for enrichment assumes the support (element) is
% divided into (k+1) sections by k cracks. It is already a general 
%implementation. It has no problem with two arbitrary fully extended
% cracks within the element. It also works for multiple branches as long as
% each branch is enriched with the correct junction function. See (Daux 2000)
r1=2*nnodes;    % uenrdofs for 1 enrichitem
r2=nnodes;      % penrdofs for 1 enrichitem
r1t=2*k*nnodes; % total uenrdofs
r2t=k*nnodes;   % total penrdofs
% 01/21/20201 Attemp to resolve the blocked pressure by existing crack. #19
if k>1
   kmat= 1000; % A big number to allow the pressure connection between cracks.
end
%% May only need one final JacobianMat, no need for discrete crack. 092920.
JMat=obj.JacobianMat; % elemental comprehensive JacobianMat
if size(JMat.Kusue,2)~=r1t % stays constant over the simulation if no new enrichitem in this element                            
    % preallocate the displacement vector for crack opening calculation
    % 02/01/2019
    JMat.Un2iEnr=zeros(r1t,1);
    % Preallocate a zero matrix for the stiffness matrix;
    Musus=zeros(r1,r1);                         % Dynamic term
    Musue=zeros(r1,r1t);
    Mueus=zeros(r1t,r1);
    Mueue=zeros(r1t,r1t);
    Kusus=zeros(r1,r1);                         % Classic stiffness matrix
    Kusue=zeros(r1,r1t);
    Kueus=zeros(r1t,r1);
    Kueue=zeros(r1t,r1t);
    Qusps=zeros(r1,r2);                         % coupling matrix
    Quspe=zeros(r1,r2t);
    Queps=zeros(r1t,r2);
    Quepe=zeros(r1t,r2t);
    Hpsps=zeros(r2,r2);                         % permeability matrix
    Hpspe=zeros(r2,r2t);  
    Hpeps=zeros(r2t,r2);  
    Hpepe=zeros(r2t,r2t);  
    Spsps=zeros(r2,r2);                         % Compressibility matrix
    Spspe=zeros(r2,r2t);
    Speps=zeros(r2t,r2);  
    Spepe=zeros(r2t,r2t);
    m=[1;1;0;1];                                % Kroneck delta in vector form
    %% Need to obtain Nuenr, Npenr, DNpenr, etc for every involved crack. 092920
    % Change: the element with multiple cracks share the same gauss points,
    % and the shape functions and their derivatives are already made
    % comprehensive by the EnrichPack.EnrichFun. 10/09/20
    GaussPnt=obj.EnrichGauss; % only one dict 10/09/20
    numgauss=length(GaussPnt);
    %% Area Integral over the element
    for igauss=1:numgauss
        Nu=GaussPnt(igauss).Nu;
        Nuenr=GaussPnt(igauss).Nuenr;
        Np=GaussPnt(igauss).Np;
        Npenr=GaussPnt(igauss).Npenr;
        DNp=GaussPnt(igauss).DNp;
        DNpenr=GaussPnt(igauss).DNpenr;
        B=GaussPnt(igauss).Bmat;
        Benr=GaussPnt(igauss).Bmatenr;
        % subarea is already included in H using subdomain method. 03302019
        dJ=1; 
%         dJ=GaussPnt(igauss).DetJ;
        H=GaussPnt(igauss).H;
        %% Compute the dynamic matrix M: Nu'*rho*Nu
        rho=Density;
        Musus=Musus+H*dJ*Nu'*rho*Nu;
        Musue=Musue+H*dJ*Nu'*rho*Nuenr;
        Mueus=Mueus+H*dJ*Nuenr'*rho*Nu;
        Mueue=Mueue+H*dJ*Nuenr'*rho*Nuenr;
        %% Compute stiffness matrix for linear elastic: B'*D*B
        D=Delastic;
        Kusus=Kusus+H*dJ*B'*D*B;
        Kusue=Kusue+H*dJ*B'*D*Benr;
        Kueus=Kueus+H*dJ*Benr'*D*B;
        Kueue=Kueue+H*dJ*Benr'*D*Benr;
        %% Compute coupling matrix Q: B'*alpha*m*Np
        alpha=Biot_alpha;
        Qusps=Qusps+H*alpha*dJ*B'*m*Np;
        Quspe=Quspe+H*alpha*dJ*B'*m*Npenr;
        Queps=Queps+H*alpha*dJ*Benr'*m*Np;
        Quepe=Quepe+H*alpha*dJ*Benr'*m*Npenr;
        %% Compute coupling matrix H: DNp'*kmat*DNp
        Hpsps=Hpsps+H*dJ*DNp'*kmat*DNp;
        Hpspe=Hpspe+H*dJ*DNp'*kmat*DNpenr;
        Hpeps=Hpeps+H*dJ*DNpenr'*kmat*DNp;
        Hpepe=Hpepe+H*dJ*DNpenr'*kmat*DNpenr;
        %% Compute coupling matrix S: Np'*Cstar*Np
        Spsps=Spsps+H*dJ*Np'*Cstar*Np;
        Spspe=Spspe+H*dJ*Np'*Cstar*Npenr;
        Speps=Speps+H*dJ*Npenr'*Cstar*Np;
        Spepe=Spepe+H*dJ*Npenr'*Cstar*Npenr;
    end
    % store the constant values to the JacobianMat(id)
    JMat.Locarray=obj.Locarray;
    JMat.LocarrayU=obj.LocarrayU;
    JMat.LocarrayP=obj.LocarrayP;
    JMat.Musus=Musus;
    JMat.Musue=Musue;
    JMat.Mueus=Mueus;
    JMat.Mueue=Mueue;
    JMat.Kusus=Kusus;
    JMat.Kusue=Kusue;
    JMat.Kueus=Kueus;
    JMat.Kueue=Kueue;
    JMat.Qusps=Qusps;
    JMat.Quspe=Quspe;
    JMat.Queps=Queps;
    JMat.Quepe=Quepe;
    JMat.Hpsps=Hpsps;
    JMat.Hpspe=Hpspe;
    JMat.Hpeps=Hpeps;
    JMat.Hpepe=Hpepe;
    JMat.Spsps=Spsps;
    JMat.Spspe=Spspe;
    JMat.Speps=Speps;
    JMat.Spepe=Spepe;
end
%% Line integral along the crack, changes every increment
% It seems better to store the line integral for each encrackitem
% separately if the leakoff calculation for each crackitem is needed.
% 10/01/20
% The elemental value remains
JMat.Kc=zeros(r1t,r1t);
JMat.Qintueps=zeros(r1t,r2);
JMat.Qintuepe=zeros(r1t,r2t);
JMat.Hintpsps=zeros(r2,r2);
JMat.Hintpspe=zeros(r2,r2t);
JMat.Hintpeps=zeros(r2t,r2);
JMat.Hintpepe=zeros(r2t,r2t);
JMat.Sintpsps=zeros(r2,r2);
JMat.Sintpspe=zeros(r2,r2t);
JMat.Sintpeps=zeros(r2t,r2);
JMat.Sintpepe=zeros(r2t,r2t);
% Loop over every Enrichitem in this element
for ienr=1:obj.EnrichNum
    % coehsive force
    Kc=zeros(r1,r1);
    % fluid pressure to the crack
    % The naming is to consistent with Queps, although Nu(std) is used
    % The enrichitem specific fractions are initiated with zeros
    Qintueps=zeros(r1,r2);                      
    Qintuepe=zeros(r1,r2t);                      
    % Leak-off terms
    Hintpsps=zeros(r2,r2);
    Hintpspe=zeros(r2,r2t);
    Hintpeps=zeros(r2t,r2);
    Hintpepe=zeros(r2t,r2t);
    Sintpsps=zeros(r2,r2);
    Sintpspe=zeros(r2,r2t);
    Sintpeps=zeros(r2t,r2);
    Sintpepe=zeros(r2t,r2t);
    linegauss=obj.LineGaussDict{ienr};
    numl=length(linegauss);
    kfi=1/gbinp.Kf;
    for igauss=1:numl
        % Retrieve the tangent and normal vector
        ds=linegauss(igauss).Ds;
        ntaud=linegauss(igauss).Ntaud*ds;
        Mtaud=linegauss(igauss).Mtaud;
        Md=Mtaud*Mtaud';
        % Preparation, the Npenr, DNpenr, should contain all enrichitems.
        % although the linegauss points are only on the current crack. The
        % number of cracks determines the repetitions of line integrals.
        % 10/01/20.
        Nu=linegauss(igauss).Nu;
%         Nuenr=linegauss(igauss).Nuenr;    
%         % Nuenr and Benr are not used in the current version 03282019
        Np=linegauss(igauss).Np;
        Npenr=linegauss(igauss).Npenr;
        DNp=linegauss(igauss).DNp;
        DNpenr=linegauss(igauss).DNpenr;
        H=linegauss(igauss).H;
        wd=linegauss(igauss).CrackOpening;
%         inigap=linegauss(igauss).MinCrackOpening;
        % note that the current crack permeability is lower than that of
        % the porous medium, which should be the reason why we observe
        % opposite direction of leak-off when injecting fluid to the
        % crack.(04242019)
        kd=(wd)^2/12/muf/1;                        % f=1.24 in [1.04,1.65]
        % For the comparison with Abaqus, f=1, and wd has an initial value
        % of 2e-4 m. 10282019
        linegauss(igauss)=linegauss(igauss).matct;
        Tc=linegauss(igauss).Tangent_coh;           % Should be obtained from matct previously.
        % start computing the matrices
        Kc=Kc+H*Nu'*Tc*Nu*ds;
        % or H*Nu'*Ntaud*Np*ds
        % to ensure sign consistency from Jacobian matrix to the Ifstd_enriched
        % Leading to symmetric LHS
%         Qintueps=Qintueps+H*skin*Nu'*ntaud*Np;
%         Qintuepe=Qintuepe+H*skin*Nu'*ntaud*Npenr;
%       leading to unsymmetric LHS
        Qintueps=Qintueps+H*Nu'*ntaud*Np;
        Qintuepe=Qintuepe+H*Nu'*ntaud*Npenr;
        %% Changed the sign of matrices related to qint on 0402 from - to + and then back to -
        % It is confirmed that qint has such a sign convention: positive
        % for flow into the domain from the crack; negative for flow
        % from the domain to the crack.(04232019)
        % NEEDS MORE ATTENTION IF DS IS NECESSARY FOR SUCH INTEGRAL WHEN MD
        % EXISTS, Changed formulation to use unit vector on 11/12/2018
        Hintpsps=Hintpsps-H*DNp'*kd*Md*DNp*wd*ds;
        Hintpspe=Hintpspe-H*DNp'*kd*Md*DNpenr*wd*ds;
        Hintpeps=Hintpeps-H*DNpenr'*kd*Md*DNp*wd*ds;
        Hintpepe=Hintpepe-H*DNpenr'*kd*Md*DNpenr*wd*ds;
        Sintpsps=Sintpsps-H*Np'*kfi*Np*wd*ds;
        Sintpspe=Sintpspe-H*Np'*kfi*Npenr*wd*ds;
        Sintpeps=Sintpeps-H*Npenr'*kfi*Np*wd*ds;
        Sintpepe=Sintpepe-H*Npenr'*kfi*Npenr*wd*ds;
    end
    % Return the normal and tangent vector to the elem-obj because 
    % Linegauss is a descendent of value class 
    obj.LineGaussDict{ienr}=linegauss;
    % Store the enrichitem specific fractions
    % Noticed it is unnecessary to store these terms after issue #30
    % because we will use elemental terms and elemental unknowns 1/1/2021
%     obj.JacobianMatDict(ienr).Kc=Kc;
%     obj.JacobianMatDict(ienr).Qintueps=Qintueps;
%     obj.JacobianMatDict(ienr).Qintuepe=Qintuepe;
%     obj.JacobianMatDict(ienr).Hintpsps=Hintpsps;
%     obj.JacobianMatDict(ienr).Hintpspe=Hintpspe;
%     obj.JacobianMatDict(ienr).Hintpeps=Hintpeps;
%     obj.JacobianMatDict(ienr).Hintpepe=Hintpepe;
%     obj.JacobianMatDict(ienr).Sintpsps=Sintpsps;
%     obj.JacobianMatDict(ienr).Sintpspe=Sintpspe;
%     obj.JacobianMatDict(ienr).Sintpeps=Sintpeps;
%     obj.JacobianMatDict(ienr).Sintpepe=Sintpepe;
    % store the elemental comprehensive matrices
    % Fractioned Kc lies along the main diagnal of elemental comprehensive Kc 
    istart=1+(ienr-1)*r1;
    iend=ienr*r1;
    JMat.Kc(istart:iend,istart:iend)=Kc;
    JMat.Qintueps(istart:iend,:)=Qintueps;
    JMat.Qintuepe(istart:iend,:)=Qintuepe;
    JMat.Hintpsps=JMat.Hintpsps+Hintpsps;
    JMat.Hintpspe=JMat.Hintpspe+Hintpspe;
    JMat.Hintpeps=JMat.Hintpeps+Hintpeps;
    JMat.Hintpepe=JMat.Hintpepe+Hintpepe;
    JMat.Sintpsps=JMat.Sintpsps+Sintpsps;
    JMat.Sintpspe=JMat.Sintpspe+Sintpspe;
    JMat.Sintpeps=JMat.Sintpeps+Sintpeps;
    JMat.Sintpepe=JMat.Sintpepe+Sintpepe;
end
%% assemble into the comprehensive Jacobian Matrix
a0=newmark.a0;
a1=newmark.a1;
a1p=newmark.a1p;
%% Prepare to rearrange the elements so that the whole Enriched Jacobian matrix WJ: 
% [Jusus,Jusue,Jusps,Juspe;
%  Jueus,Jueue,Jueps,Juepe;
%  Jpsus,Jpeue,Jpsps,Jpspe;
%  Jpeus,Jpeue,Jpeps,Jpepe];
% The actual WJ in the end is reordered to be consistent with the dofs [us;ps;ue;pe]
Jusus=-a1*(a0*JMat.Musus+JMat.Kusus);
Jusue=-a1*(a0*JMat.Musue+JMat.Kusue);
Jusps=a1*JMat.Qusps;
Juspe=a1*JMat.Quspe;
Jueus=-a1*(a0*JMat.Mueus+JMat.Kueus);
Jueue=-a1*(a0*JMat.Mueue+JMat.Kueue+JMat.Kc);
% when it comes to a new increment the Kc should be left out because the
% unkonwns is total displacement but Kc is based on incremental
% displacement. The cohesion traction(from the previous increment) was
% already added in the intforcer. 04122019
Jueuenew=-a1*(a0*JMat.Mueue+JMat.Kueue);    
Jueps=a1*(JMat.Queps+JMat.Qintueps);
Juepe=a1*(JMat.Quepe+JMat.Qintuepe);
Jpsus=a1*JMat.Qusps';
Jpsue=a1*transpose(JMat.Queps+JMat.Qintueps);
% IMPORTANT ERROR: SIGN OF HINTPSPS, SINTPSPS. (FIXED ON 04022019)
Jpsps=JMat.Hpsps-JMat.Hintpsps+a1p*(JMat.Spsps-JMat.Sintpsps);
% IMPORTANT ERROR: SIGN OF HINTPSPE, SINTPSPE. (FIXED ON 04022019)
Jpspe=JMat.Hpspe-JMat.Hintpspe+a1p*(JMat.Spspe-JMat.Sintpspe);
Jpeus=a1*JMat.Quspe';
Jpeue=a1*transpose(JMat.Quepe+JMat.Qintuepe);
% IMPORTANT ERROR: SIGN OF HINTPEPS, SINTPEPS. (FIXED ON 04022019)
Jpeps=JMat.Hpeps-JMat.Hintpeps+a1p*(JMat.Speps-JMat.Sintpeps);
% IMPORTANT ERROR: SIGN OF HINTPEPE, SINTPEPE. (FIXED ON 04022019)
Jpepe=JMat.Hpepe-JMat.Hintpepe+a1p*(JMat.Spepe-JMat.Sintpepe);
% The whole elemental Jacobian Matrix
JMat.WJ=[Jusus,Jusps,Jusue,Juspe;
        Jpsus,Jpsps,Jpsue,Jpspe;
        Jueus,Jueps,Jueue,Juepe;
        Jpeus,Jpeps,Jpeue,Jpepe];
JMat.WJnew=[Jusus,Jusps,Jusue,Juspe;
        Jpsus,Jpsps,Jpsue,Jpspe;
        Jueus,Jueps,Jueuenew,Juepe;
        Jpeus,Jpeps,Jpeue,Jpepe];
end

