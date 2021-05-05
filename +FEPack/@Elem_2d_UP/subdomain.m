function subdomain(obj,varargin)
%SUBDOMAIN use the LocalInt to create triangular subdomains for integral
%purpose. This method should be called by higher level class when every
% EnrCrackBody has finished the enriching. Otherwise, it is just a waste
% of computation cost.10/31/2020

% Explanation of outputs:
% 1. gp: the xi and eta in terms local parent element coordinates
% 2. gw: the weights of the integration points
% Explanation of the inputs:
% 1. obj is the Elem_2d_UP element
% 2. the number of gauss points within one triangle (p<=13)

%   The seeds in terms of local coordinates, xi and eta have been created
%   and assigned from ToolPact.OpenGeo.intersection method 03302019.

%% 10/30/20 Generate the polygons for subdomains
% Now the seeds will be generated here using new method for more
% general cases when there are more than one crcks inside this element.
% This method should be better than the previous hand-coded criterion to
% find sub-polygons. Note convex hull method is not necessary in this
% context as we would have the vertices on the boundaries but not internal
% point cloud.

% 1.  Read the cracks intersections with the boundaries (obj.LocalInt) and
% express the crack segements in functions
NLines=obj.EnrichNum;
RealEnrichLines=find(~obj.Smeared); % should not be blocked in the commit
% d5fe520e92eda7eaa5f624205294850159b9196b
MaxPart=NLines*(NLines+1)/2+1; % Maximum possible partitions from NLines
polygons=cell(1,MaxPart);
localints=zeros(NLines*2,2);
localints_inrow=zeros(NLines,4);
LineHandles=cell(1,NLines); % a cell array to store the function handles
% create line function handles by two points
for i=1:NLines
    % BUG after introducing smeared and realenrich, fixed on 03/16/21.
%     can be done alternatively iline=obj.get_realenrichind(i);
    iline=RealEnrichLines(i); % TO find the not smeared line index.
    x1=obj.LocalInt{iline}(1,1);
    y1=obj.LocalInt{iline}(1,2);
    x2=obj.LocalInt{iline}(2,1);
    y2=obj.LocalInt{iline}(2,2);
    localints(2*i-1:2*i,:)=obj.LocalInt{iline};
    localints_inrow(i,:)=[x1,y1,x2,y2];
    LineHandles{i}=@(x,y) y-y1-(y2-y1)*(x-x1)/(x2-x1);
end

% 2. Find the intersections of these crack segments and store the
% intersections together with the localpnts, nodes as the vertices for the
% sub-polygons
intersections=[];
if NLines>1
    AllPairs=sortrows(combnk(1:NLines,2));
    % LineIntersection may not work because it does not constrain it for this
    % segment.
    results = lineSegmentIntersect(localints_inrow(AllPairs(:,1),:) , localints_inrow(AllPairs(:,2),:) );
    logicalind=diag(results.intAdjacencyMatrix);
    px=diag(results.intMatrixX);
    px=px(logicalind);
    py=diag(results.intMatrixY);
    py=py(logicalind);
    intersections=[px,py]; % intersections of the crack segments
end
nodes=[-1,-1;1,-1;1,1;-1,1];
% 11/23/20 The unique may be too strong because some points may be very
% close to another point. Change it to uniquetol for desired result.
allvertices=uniquetol([nodes;localints;intersections],1e-5,'ByRows',true);

% 3. Plug the vertices into the equations and store the sign in a row
% vector of the size crack segments. >=0: 1, <=0 : 0. stored in signs=[1 0]
% real sign determined by sign(linehandle)
signmat=zeros(size(allvertices,1),NLines);
% group mat will furtehr determine if the signmat is >=0 and if it is <=0 as well
groupmat=true(size(signmat,1),2,NLines);
tol=1e-6; % The numerical calculation of the intersection point may be off the lines
for iline=1:NLines
%     iline=RealEnrichLines(i); % TO find the not smeared line index.
    temp=LineHandles{iline}(allvertices(:,1),allvertices(:,2));
    temp(abs(temp)<tol)=0;
    signmat(:,iline)=sign(temp);                % real sign matrix
    groupmat(:,1,iline)=signmat(:,iline)>=0;    % include bounary
    groupmat(:,2,iline)=signmat(:,iline)<=0;
end

% 4. Group the vertices by groupmat pair. Then we would have
% the vertices for each subdomain (polygon).
npart=0;
% 11/22/20 Fix on the single line mode
if NLines == 1    % There will be only two polygons created by one line
    npart=2;
    polygons{1}=find(groupmat(:,1));
    polygons{2}=find(groupmat(:,2));
elseif NLines > 1 % find the polygons by looping through AllPairs
    for ipair=1:size(AllPairs,1)
        line1=AllPairs(ipair,1);
        line2=AllPairs(ipair,2);
        for i=1:2
            for j=1:2
                group=[groupmat(:,i,line1),groupmat(:,j,line2)];
                part=find(all(group,2)); % treat all logical on each row
                if length(part)>2 % two more points to make a closed 2d shape
                    npart=npart+1;
                    polygons{npart}=part;
                end
            end
        end
    end
end
% Check if the subdomain is already created (by the other crack)
% This is very important, otherwise, the later call of this method
% will erase the previous gaussian points with enrichment info. 11/23/2020
if npart~=obj.PolygonNum
    obj.PolygonNum=npart;
    polygons=polygons(1:npart); % () bracing access the sets of cell arrays, {}access the contents
    iniseeds=cell(1,npart);
    % 5. Based on the area of subdomains, choose how many inner points will be
    % generated for the delaunay triangulation. Store seeds in cell arrays.
    for ipart=1:npart
        tempseeds=[allvertices(polygons{ipart},1),allvertices(polygons{ipart},2)]; %(xs,ys)
        % do not use polyarea here as it requires the vertices at the circular order
        % use convhull to get the convex hull of the vertices in the right
        % order
        [k,area]=convhull(tempseeds(:,1),tempseeds(:,2));
        tempseeds=tempseeds(k(1:end-1),:); % this is very important
        if area>2
            vertices_x=tempseeds(:,1);
            vertices_y=tempseeds(:,2);
            tri=delaunay(vertices_x,vertices_y);
            tri_x=vertices_x(tri);
            tri_y=vertices_y(tri);
            centroid=[sum(tri_x,2)/3,sum(tri_y,2)/3];
        elseif area>0.25
            centroid = ploygon_centroid(tempseeds(:,1),tempseeds(:,2)); % get the centroid
        else % do not need further triangulation
            centroid=[];
        end
        iniseeds{ipart}=[tempseeds;centroid];
    end
    
    %% Generate the gaussian points for each subdomain
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
    % id=obj.Enrich==id;
    % seeds_all=obj.Seeds;        % all subpolygons for current crack id
    GNcoord=[obj.X,obj.Y];      % global coordinates of the nodes of the element
    GP=zeros(100,2); % preallocate
    GW=zeros(100,1);
    tpt=0;
    for ipart=1:npart
        seeds=iniseeds{ipart};
        % recover the global coordinates GScoord of the seeds
        Nquad= @(xi,eta) [(1-xi).*(1-eta),(1+xi).*(1-eta),(1+xi).*(1+eta),(1-xi).*(1+eta)]/4;
        Nquadmat=Nquad(seeds(:,1),seeds(:,2));
        GScoord=Nquadmat*GNcoord;
        % triangulation
        warning('off', 'MATLAB:delaunay:DupPtsDelarunayWarnId')
        tri = delaunay(seeds(:,1),seeds(:,2));
%             triplot(tri,seeds(:,1),seeds(:,2));
%             triplot(tri,GScoord(:,1),GScoord(:,2));
%             hold on
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
				% Checked with Serfeim_bakalakos thesis 2000.
                subarea=polyarea(gcoord(:,1),gcoord(:,2)); %AK
                % polyarea used here has no problem because triangular
                % vertices can be arbitrarily ordered.
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
end

function centroid = ploygon_centroid(x,y)
xs = circshift(x,-1);
ys = circshift(y,-1);
area = 0.5*sum (x.*ys-xs.*y); % also require the right ordering
xc = sum((x.*ys-xs.*y).*(x+xs))/(6*area);
yc = sum((x.*ys-xs.*y).*(y+ys))/(6*area);
centroid=[xc,yc];
end

