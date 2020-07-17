function crtstif( obj, newmark, iinc, gbinp, blending )
%CRTST method of elem_2d_EP to calculate the stiffness matrix using gauss
%point class
%   This function is slightly different from the calstifmatrix in class
%   elem_2d due to the new member GaussPntDict in the elem_2d_EP
%   Each gauss point object calculates the consistent tangent operator,D
%   matrix. This function is more like a assembler.
% global  Delastic Density Biot_alpha kmat Cstar; 
if blending
    D=gbinp.Delastic;
    rho=gbinp.Density;
    alpha=gbinp.Biot_alpha;
    kmat=gbinp.kmat_crack;
    Cstar=gbinp.Cstar;
else
    D=gbinp.Delastic;
    rho=gbinp.Density;
    alpha=gbinp.Biot_alpha;
    kmat=gbinp.kmat;
    Cstar=gbinp.Cstar;
end
if blending || iinc==1                            % compute all components for blending elements
    nnodes=obj.NoNodes;
    r1=2*nnodes;
    r2=nnodes;
    % Preallocate a zero matrix for the stiffness matrix;
    Melem=zeros(r1,r1);                         % Dynamic term
    Kelem=zeros(r1,r1);                         % Classic stiffness matrix
    Qelem=zeros(r1,r2);                         % coupling matrix
    Helem=zeros(r2,r2);                         % permeability matrix
    Selem=zeros(r2,r2);                         % Compressibility matrix
    m=[1;1;0;1];                                % Kroneck delta in vector form
    numgauss=length(obj.GaussPntDictM);
    for igauss=1:numgauss
        % calculate the tangent operator consistent with the last converged
        % solution: modified newton-raphson
        Nu=obj.GaussPntDictM(igauss).Nu;
        Np=obj.GaussPntDictM(igauss).Np;
        DNp=obj.GaussPntDictM(igauss).DNp;
        B=obj.GaussPntDictM(igauss).Bmat;
        dJ=obj.GaussPntDictM(igauss).DetJ;
        H=obj.GaussPntDictM(igauss).H;
        %% Compute the dynamic matrix M: Nu'*rho*Nu
        if obj.NoNodes==3 || obj.NoNodes==6
            Melem=Melem+0.5*H*dJ*Nu'*rho*Nu;            % The formula of area integral over a triangle has this coefficient 1/2.
        else
            Melem=Melem+H*dJ*Nu'*rho*Nu;
        end
        %% Compute stiffness matrix for linear elastic: B'*D*B
        if obj.NoNodes==3 || obj.NoNodes==6
            Kelem=Kelem+0.5*H*dJ*B'*D*B;            % The formula of area integral over a triangle has this coefficient 1/2.
        else
            Kelem=Kelem+H*dJ*B'*D*B;
        end
        %% Compute coupling matrix Q: B'*alpha*m*Np
        if obj.NoNodes==3 || obj.NoNodes==6
            Qelem=Qelem+0.5*H*alpha*dJ*B'*m*Np;            % The formula of area integral over a triangle has this coefficient 1/2.
        else
            Qelem=Qelem+H*alpha*dJ*B'*m*Np;
        end
        %% Compute coupling matrix H: DNp'*kmat*DNp
        if obj.NoNodes==3 || obj.NoNodes==6
            Helem=Helem+0.5*H*dJ*DNp'*kmat*DNp;            % The formula of area integral over a triangle has this coefficient 1/2.
        else
            Helem=Helem+H*dJ*DNp'*kmat*DNp;
        end
        %% Compute coupling matrix S: Np'*Cstar*Np
        if obj.NoNodes==3 || obj.NoNodes==6
            Selem=Selem+0.5*H*dJ*Np'*Cstar*Np;              % The formula of area integral over a triangle has this coefficient 1/2.
        else
            Selem=Selem+H*dJ*Np'*Cstar*Np;
        end
    end
    obj.M=Melem;
    obj.K=Kelem;
    obj.Q=Qelem;
    obj.H=Helem;
    obj.S=Selem;
else                                                        % For the subsequent increment, only read components from object
    Melem=obj.M;
    Kelem=obj.K;
    Qelem=obj.Q;
    Helem=obj.H;
    Selem=obj.S;
end
%% assemble into the comprehensive Jacobian Matrix
a0=newmark.a0;
a1=newmark.a1;
a1p=newmark.a1p;
obj.StifMatrix=[-a0*a1*Melem-a1*Kelem, a1*Qelem; a1*Qelem', Helem+a1p*Selem];
% %% Prepare to rearrange the elements so that the whole StifMatrix will be in the order of [ux,uy,p...]
% obj.JM11=-a0*a1*Melem-a1*Kelem;   
% obj.JM12=a1*Qelem;
% obj.JM21=a1*Qelem';
% obj.JM22=Helem+a1p*Selem;
end

