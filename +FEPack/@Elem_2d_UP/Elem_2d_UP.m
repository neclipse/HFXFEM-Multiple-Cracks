classdef Elem_2d_UP < handle
    properties
        dXn2i                   % (standard) unknowns solved from current iteration [dx1,dy1,dp1,dx2,dy2,dp2...]	
        Un2i                    % (standard) Un1+dUn2i
        Pn2i                    % (standard) Pn1+dPu2i
        Un2t1
        Un2t2
        Pn2t1
        LocalInt=cell(1,3);   % local coordinates of the intersections   
        GlobalInt=cell(1,3);  % global coordinates of the intersections
        LineGaussDict=cell(1,3);      % object array of the GaussPnt for line integral only
        EnrichGauss;                  % The EnrichGauss points with assembled Nuenr, Npenr, and their derivatives
    end
    properties (SetAccess = private, GetAccess = public)
        Ind
        Type                    % Element type reheritated from the class Domain
        MatType
        NodList                 % Global order of the nodes in the element
        NodDict
        X
        Y
        Area
        GaussPntDictM;          % object array of class standard GaussPnt
        NoNodes                 % Number of total nodes in this element
        IntLoadVec              % Element internal load vector
        ExtLoadVec              % Element external load vector
        StifMatrix              % Coupled stiffness matrix for the coupled equation
        M
        K
        Q
        H
        S
        Length             		% length of each side
        Locarray                % global location array of mixed dofs,[ux,uy,p]
        LocarrayU               % global location array of only U, [ux,uy],  not obsolete, used in linsys.upconf
        LocarrayP               % global location array of only p, [p]
        EnrichNum=0;           % Number of Enrichment items involved with this element
        Enrich=zeros(1,3);    % Enrichment flags: special structure:[id1,id2,id3...]
        Smeared=false(1,3);   % Smeared flags: corresponding to the obj.Enrich.
        PolygonNum=0;         % number of sub polygons created by the cracks. default=1, to keep track of subdomain
        JacobianMat           % The comprehensive elemental JocobianMat, object of FEPack.JacobianMat
        JacobianMatDict       % Enriched Jacobian matrix for the enriched element, objects of FEPack.JacobianMat
        Stress                 % Total stresses at the element centroid
        Stressp                % Effective stresses at the element centroid
    end
    properties (Dependent)
        RealEnrich             % obj.Enrich(~obj.Smeared)
    end
    
    methods
        function obj=Elem_2d_UP(index,elemtype,mattype,nodlist,noddict,gaussdictm)
            if nargin>0
               obj.Ind=index;
               obj.Type=elemtype;
               obj.MatType=mattype;
               obj.NodList=nodlist; 
               % nodlist is ordered in counterclockwise by vertices first and side nodes last
               obj.NodDict=noddict;  
               obj.X=[noddict.X]';
               obj.Y=[noddict.Y]';
               obj.NoNodes=length(noddict);
%                obj.GaussPntDictM=gaussdictm;
               for igauss=1:length(gaussdictm)
                  gaussdictm(igauss)=gaussdictm(igauss).preparing(obj.X,obj.Y); 
               end
               obj.GaussPntDictM=gaussdictm;
            end
        end
		function dislist = callength(obj) 
            if obj.NoNodes == 3||4
                xlist1=obj.X;
                xlist2=[obj.X(2:end);obj.X(1)];
                ylist1=obj.Y;
                ylist2=[obj.Y(2:end);obj.Y(1)];
                xdis=xlist2-xlist1;
                ydis=ylist2-ylist1;
                dislist=sqrt(xdis.^2+ydis.^2);
                obj.Length=dislist;
            else
                error('not supported yet');
            end
        end
        
        function setenrich(obj,id,initialmode)
            % obj.Enrich stores all the ids of involved enrichitem with
            % this element
            L=id~=obj.Enrich;
            if isempty(obj.JacobianMatDict)
                % preallocate the Jacobian matrices for potential cracks
                % CHANGE SEMICOLON TO COMMA AS SEMICOLON WOULD GIVE A DICTIONARY OF ALL
                % REPEATED HANDLE 11/16/2018.
                JMatDict(1,3)=FEPack.JacobianMat();
                obj.JacobianMatDict=JMatDict;
                % prepare the assembled elemental Jacobian matrix
                obj.JacobianMat=FEPack.JacobianMat();
            end
            if all(L)
                % the indexing is not good 092820, to change
                % Changed on 10/02/20
                newnum=nnz(obj.Enrich)+1;
                if initialmode~=2 % not "smeared crack mode"
                    obj.EnrichNum=obj.EnrichNum+1;
                    obj.Enrich(newnum)=id;
                    obj.JacobianMatDict(newnum)=FEPack.JacobianMat(id);
                elseif initialmode==2  %'smeared crack'
                    obj.Enrich(newnum)=id;
                    obj.Smeared(newnum)=true;
                    obj.JacobianMatDict(newnum)=FEPack.JacobianMat(id);
                end
            end
        end
        function real_ind=get_realenrichind(obj,i)
            % i is the numerical index for obj.Enrich
            % real_ind is the correct index for the real enrcrack in obj.Enrich 
            notsmeared=find(~obj.Smeared);
            real_ind=notsmeared(i);
        end
        function value=get.RealEnrich(obj)
%             value=zeros(size(obj.Enrich));
%             value(~obj.Smeared)=obj.Enrich(~obj.Smeared);
              value=obj.Enrich(~obj.Smeared); % Enrich=[2,1,0], smeared=[1,0,0], realenrich=[1,0]
        end
        function opensmeared(obj,enrichid)
            % update the flags and EnrichNum when the element opens in
            % enrcrack.postprocess
            obj.EnrichNum=obj.EnrichNum+1;
            obj.Smeared(obj.Enrich==enrichid)=false;
        end
        %% Function Prototyping
        assemble_locarray( obj );                           % assemble enriched locarray from JacobianMatDict to obj.JacobianMat.
        crtstif(obj,newmark,iinc,gbinp, blending);           % create element stiffness matrix, call matct in gausspnt
        crtstif_enriched_elemental( obj, newmark, gbinp, blending);
%         crtstif_enriched_1Dflow( obj, newmark, id , gbinp); % 1D crack fluid flow adopted and compressibility of fracking fluid ignored
        givelocarray(obj,varargin);
        givelocarray_enriched(obj,crackid,varargin);
		load=ifstd(obj,newmark);     			% compute internal force vectors; call listra and matsu
        [IntLoadAll,stagechangeflag]=ifstd_enriched( obj,newmark, stagecheck);
        calarea(obj);
		[flagie,flagi,flage,flagoe,area] = isinside_vec(obj,plist);
		subdomain(obj, varargins);
        linegauss( obj,id,cohesive,perforated,varargin );
        [stressp, stress]=extraplstress( obj, xi, eta);
        calstress(obj,storage);
        calstress_enriched(obj,storage);
    end
    
end