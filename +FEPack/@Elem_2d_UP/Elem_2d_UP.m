classdef Elem_2d_UP < handle
    properties
        Ind
        Type                    % Element type reheritated from the class Domain
        MatType
        NodList                 % Global order of the nodes in the element
        NodDict
        GaussPntDictM;                 % object array of class standard GaussPnt
        LineGaussDict=cell(1,3);      % object array of the GaussPnt for line integral only
        EnrichGauss;                  % The EnrichGauss points with assembled Nuenr, Npenr, and their derivatives
        EnrichGaussDict=cell(1,3);    % object array of GaussPnt for enriched element 2d integral, specific to each enritem
        X
        Y
        Area
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
        dXn2i                   % (standard) unknowns solved from current iteration [dx1,dy1,dp1,dx2,dy2,dp2...]	
        Un2i                    % (standard) Un1+dUn2i
        Pn2i                    % (standard) Pn1+dPu2i
        Un2t1
        Un2t2
        Pn2t1
        Stress                 % Total stresses at the element centroid
        Stressp                % Effective stresses at the element centroid
        EnrichNum=0;          % The total involved enriched items in this element
        Enrich=zeros(1,3);    % Enrichment flags: special structure:[id1,id2,id3...]
        Seeds;                % seeds points to create triangular subdomain, only store the latest updated info.
        LocalInt=cell(1,3);   % local coordinates of the intersections   
        GlobalInt=cell(1,3);  % global coordinates of the intersections
        JacobianMat           % The comprehensive elemental JocobianMat, object of FEPack.JacobianMat
        JacobianMatDict       % Enriched Jacobian matrix for the enriched element, objects of FEPack.JacobianMat
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
        
        function setenrich(obj,id)
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
               obj.EnrichNum=obj.EnrichNum+1; 
               obj.Enrich(obj.EnrichNum)=id; 
               obj.JacobianMatDict(obj.EnrichNum)=FEPack.JacobianMat(id);
               
%                NOT necessary after the revision on 02262019, stdnodes are
%                now stored in the "mygeo" of the encrack.
%                stdnodes=1:length(obj.NodList);
%                L=false(1,length(obj.NodList));
%                for iN=1:length(obj.NodList)
%                    if all(obj.NodDict(iN).Enrich~=id)
%                        L(iN)=true;
%                    end
%                end
%                stdnodeslist=stdnodes(L);
%                obj.NodstdList{id}=stdnodeslist;
           end
        end
        
        %% Function Prototyping
        assemble_locarray( obj );                           % assemble enriched locarray from JacobianMatDict to obj.JacobianMat.
        crtstif(obj,newmark,iinc,gbinp, blending);           % create element stiffness matrix, call matct in gausspnt
        crtstif_enriched(obj,newmark,id,gbinp,blending);    % create element stiffness matrix, call matct in gausspnt, for enrichitem id.
        crtstif_enriched_elemental( obj, newmark, gbinp, blending);
%         crtstif_enriched_1Dflow( obj, newmark, id , gbinp); % 1D crack fluid flow adopted and compressibility of fracking fluid ignored
        givelocarray(obj,varargin);
        givelocarray_enriched(obj,crackid,varargin);
		load=ifstd(obj,newmark , calstress);     			% compute internal force vectors; call listra and matsu
        [IntLoadAll,stagechangeflag]=ifstd_enriched( obj,newmark, stagecheck , calstress);
        calarea(obj);
%         flag=isinside(obj,x,y);				% check if point(x,y) is inside or on the edge of the element
		[flagie,flagi,flage,flagoe,area] = isinside_vec(obj,plist);
		subdomain(obj, varargins);
        linegauss( obj,id,cohesive,perforated,varargin );
        [stressp, stress]=extraplstress( obj, xi, eta);
        calstress(obj);
        calstress_enriched(obj);
    end
    
end