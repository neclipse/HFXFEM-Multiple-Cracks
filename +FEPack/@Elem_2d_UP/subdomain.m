function [gp,gw]=subdomain(obj,id,varargin)
%SUBDOMAIN use the seeds to create triangular subdomains for integral
%purpose. called by EnrCrackBody. Last update date:03302019.
%   The seeds in terms of local coordinates, xi and eta have been created
%   and assigned from ToolPact.OpenGeo.intersection method
% Explanation of outputs:
% 1. gp: the xi and eta in terms local parent element coordinates
% 2. gw: the weights of the integration points
% Explanationo f the inputs:
% 1. obj is the Elem_2d_UP element
% 2. the number of gauss points within one triangle (p<=13)
if ~isempty(varargin)
    p=varargin{1};
else
    p=3;
end
% obtain the local parent coordinates of all points to create delaunay
% triangles note seeds_all is still a cell array because obj.Seeds has two
% levels: first level is the crack id and the second level is the
% subdomains for triangulation 03282019
% use logical array id to relieve the single index id. 10/02/20
id=obj.Enrich==id;
seeds_all=obj.Seeds{id};    % all subpolygons for current crack id
GNcoord=[obj.X,obj.Y];      % global coordinates of the nodes of the element
GP=zeros(60,2); % preallocate
GW=zeros(60,1);
tpt=0;
for isub=1:length(seeds_all)
    seeds=seeds_all{isub};
    % recover the global coordinates GScoord of the seeds
    Nquad= @(xi,eta) [(1-xi).*(1-eta),(1+xi).*(1-eta),(1+xi).*(1+eta),(1-xi).*(1+eta)]/4;
    Nquadmat=Nquad(seeds(:,1),seeds(:,2));
    GScoord=Nquadmat*GNcoord;
    % triangulation
    warning('off', 'MATLAB:delaunay:DupPtsDelarunayWarnId')
    tri = delaunay(seeds(:,1),seeds(:,2));
%     triplot(tri,seeds(:,1),seeds(:,2));
%     triplot(tri,GScoord(:,1),GScoord(:,2));
%     hold on
    [q,w] = gausstable(p,'TRI');
    NPT=p*size(tri,1);% number of total integration poionts for this subpolygon
    gp=zeros(NPT,2);
    gw=zeros(NPT,1);
    % Loop over subtriangles to get the quadrature points and their weights
    pt=1;
    for e=1:size(tri,1)
        % in local coordinates (Xi, Eta) for the quadrilateral element
        coord = seeds(tri(e,:),:);  
        gcoord= GScoord(tri(e,:),:);% global coordinates for the subtriangle
        for i = 1:length(w)
            % local coordinates (ksi, ita) for the triangle element
            ksi=q(i,1);              
            ita=q(i,2);
            %local coordinates (Xi, Eta) of the gauss points within the
            %subtriangle for the quadrilateral element
            Ntri=[1-ksi-ita,ksi,ita];
            gp(pt,:)=Ntri*coord;           
%             subarea=polyarea(coord(:,1),coord(:,2));   
            % the gauss quadrature over an triangular element has 0.5 infront
            % of the w(i)*f(ksi,ita) comparing to the quadrilateral element
            % because the area of the unit triangle is 0.5 (03282019)
%             gw(pt)=w(i)*subarea;  %fully expressed by (w(i)*1/2)*area/(1/2); 
            % the area of the subtriangle should be used to adjust the
            % integration weights (area of unit square)/(area of unit
            % triangle) * (area of subtriangle)/(area of unit square)
            % =(area of subtriangle)/ (area of unit triangle) --03292019
            % went back to this version after realising that the follwoing
            % discussion forgot to account for the detJ in the gauss
            % quadrature for 2d. But it is still not clear if the detJ is
            % being correctly calculated in the current form 03302019.
%             gw(pt)=w(i)/2;  

            % UPDATED version 03302019. SHOULD BE CORRECT USING THE
            % SUBAREA OF THE SUBTRIANGLE IN THE GLOBAL COORDINATES, AK. THE
            % INTEGRATION FORMULA FOR INT(F(x,y))BECOMES 
            % SUM[2*Ak*w(i)/2*F(P(xi_i,eta_i),Q(xi_i,eta_i))]. Note detJ
            % is already expressed by 2*AK. so the easy way to accomadate
            % the change is to combine the w(i)/2 and 2*Ak for gw(pt).
            subarea=polyarea(gcoord(:,1),gcoord(:,2)); %AK
            gw(pt)=w(i)*subarea;
            pt=pt+1;
        end
    end
    GP(tpt+1:tpt+NPT,:)=gp;
    GW(tpt+1:tpt+NPT)=gw;
    tpt=tpt+NPT;
end
% remove preallocated zeros
GP=GP(1:tpt,:);
GW=GW(1:tpt);
nnodes=length(obj.NodList);
% NOTE the order of the gauss points in these enriched elements
%are not consistent with those in traditional elements (1,3,4,2)
gbinp=obj.GaussPntDictM(1).GBINP;
gaussdictm(1,tpt)=FEPack.GaussPnt_LE_UP(); % generate an void object array
    for igauss=1:tpt
        gaussdictm(igauss)=FEPack.GaussPnt_LE_UP(GP(igauss,1),GP(igauss,2),GW(igauss),nnodes,gbinp);
        gaussdictm(igauss)=gaussdictm(igauss).preparing(obj.X,obj.Y,obj.EnrichNum);
    end
% replace the GaussPntDictM with the newly generated gaussdictm
% obj.EnrichGaussDict{id}=gaussdictm;
% The comprehensive gaussdict should be also updated. Actuall, the
% comprehensive gaussdict should be the same as every cell within the
% obj.EnrichGaussDict after update. 10/16/20
obj.EnrichGauss=gaussdictm; 
end

