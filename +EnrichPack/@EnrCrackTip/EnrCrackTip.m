classdef EnrCrackTip < EnrichPack.EnrichItem
    properties (Access = public)
        Itip                                                                 % index of the tip of the current crack, 1-- the last s, 2-- the first
        Xct                                                                  % vector, history of crack x-coordinates of the crack tip, the current position is the last component
        Yct                                                                  % vector, history of crack Y-coordinates of the crack tip, the current position is the last component
        Ki                                                                   % vector,stress intensity factor of Mode I
        Kii                                                                  % vector,stress intensity factor of Mode II
        Omega                                                                % vector,tilting angle of crack tip, counterclockwisely from horizontal
        Radius                                                               % parameter l for the gaussian weight w, also assumed as the radius of interaction 05172019
        Stressp                                                              % the stress tensonr (in vector form) at the crack tip, maybe used to determine the growing
        Growcheck        % An object of ToolPack.PropagationCheck
        Growdirection    % An object of ToolPack.CrackGrowDirectionLaw (Discarded in the version of Maxps, 08132019)
        Growincrement    % Grow increment calculated by various methods, here we can use obj.elementincremnt;
        NextElem         % Element object
        NewNodes         % Newly added nodes to be enriched
        NewElems         % Newly added elements to be enriched.
        UpdateElems      % May include the old tip element to be updated, the enrichement and the stiffness matrix
    end
    properties (Access = protected)
        Isactive = true
        Mycrackbody
    end
    
    methods
        
        function obj = EnrCrackTip(type,elemdict,nodedict,mygeo,itip)
            obj = obj@EnrichPack.EnrichItem(type,elemdict,nodedict);
            obj.Mygeo=mygeo;
            obj.Itip=itip;
            obj.Mesh=mygeo.Mesh;
            obj.initiate;
        end
        function [growflag,unstablegrow,cutflag]=lookahead(obj)
            % An important method for the crack tip for crack propagation:
            % (in the explicit crack front tracking scheme) the method will
            % called at the end of increment right after the
            % encrack.postprocess. For now, there is no need to recheck the
            % equilibrium as the initial enriched dofs will be zeros for the
            % new enriched elements. NOT CORRECT. If the initial cohesion is
            % not identical to the tensile traction before inserting the
            % crack, the equilibrium will be disturbed. 06172019 However,
            % the recheck of equilibrium can be postponed if only one
            % element will be inserted at most because we do not need
            % accurate stress update to determine if the crack will continue
            % to grow.
            growflag=false;
            unstablegrow=false;
            cutflag=false;
            % Isactive used to control if look ahead is really implemented
            % 12/07/2020
            if obj.Isactive
                obj.findnextelem; % only find the nextelem and find the intersections
                % but do not extend the crack tip into the nextelem
                % If later on the grow check needs to check if obj.NextElem
                % is already cut by other cracks, do the check here.
                % 11/07/20
                if obj.NextElem.EnrichNum>0
                    % 1. May use different growcheck criterion
                    fprintf('The NextElem of tip %d of crack %d has existing cracks.\n',obj.Itip, obj.Id);
                    % Temporarily increase the tolerance to 100% to allow
                    % the crack extend to the next element.(not adopted)
                    obj.Growcheck.Mode='all';
%                     obj.Growcheck.Tol2=0.5;
                    % may use nextelem.stressp determine the possible
                    % interatcion scenario.
                    % if obj.NextElem.Streesp
                else
                    % 2. obj.GrowCheck.growcheck;
%                     obj.Growcheck.Tol2=0.15;
                    obj.Growcheck.Mode='tip';
                end
                obj.Growcheck.cal_pvariable(obj.Stressp,obj.Itip, obj.Omega,obj.NextElem);
                obj.Growcheck.growcheck;  % unstable flag can be the output
                growflag=obj.Growcheck.Growflag;
                unstable=obj.Growcheck.Unstable;
                % The current time increment is too large for a good crack
                % growth prediction, return to cut back the current inc.
                % Added on 07112019
                if unstable==true
                    if growflag==false
                        cutflag=true;
                        return;
                    else
                        unstablegrow=true;
                    end
                end
            end
        end
        function checkactive(obj)
            % The Isactive is primarily determined by Mygeo.Rtip 12/07/20
            if all(obj.Mygeo.Rtips~=obj.Itip)
                obj.Isactive = false;
            end
        end
        function setinteractedelem(obj,itip)                                     % Find the interacted elements and assign to the obj.InteractedElem, implemented by mygeo
            obj.Interactedelem=obj.Mygeo.Rtipelements(itip);
            obj.INTELEM=obj.Elemdict(obj.Mygeo.Rtipelements(itip));
        end
        function setenrichednode(obj)
            obj.Enrichednode=obj.INTELEM.NodList;
            obj.ENNODE=obj.INTELEM.NodDict;
        end
        
        function [nextelem,pnts,localpnts,Phi]=findnextelem(obj,varargin)
            if obj.Isactive
                % Find next good element
                if ~isempty(varargin)
                    gtheta=varargin{1};
                    % to handle the second tip, 08162019. 
                    % not needed after psu_exact is being used, 11192019
%                     if obj.Itip==2
%                         gtheta=pi+gtheta;
%                     end
                else
                    gtheta=obj.Omega;
                end
                goodelem=false;
                initial=[obj.Xct,obj.Yct];
                currentelem=obj.INTELEM;
                while ~goodelem
                    %                elems=obj.Mygeo.findelems(initial,'in_edge');
                    %% IMPORTANT BUG: SHOULD NOT USE THIS METHOD TO FIND NEXTELEM
                    % BECAUSE IT IS NOT GUARANTEED THAT THE CRACK WILL GROW TO
                    % THE NEIGHBORING ELEMENT SHARING THE TIP EDGE. 08182019
                    % (not well fixed)
                    %                if length(elems)==2
                    %                    nextelem=elems(elems~=currentelem.Ind);
                    %                elseif length(elems)==4
                    %                    r=0.2*min(currentelem.Length);
                    % %                    gtheta=theta+obj.Omega;
                    %                    xtemp=obj.Xct+r*cos(gtheta);
                    %                    ytemp=obj.Yct+r*sin(gtheta);
                    %                    nextelem=obj.Mygeo.findelems([xtemp,ytemp],'in_edge');
                    %                end
                    ratio=0.001;
                    % Loop to find the next element along gtheta.
                    nextelem=currentelem.Ind;
                    while nextelem==currentelem.Ind
                        r=ratio*min(currentelem.Length);
                        xtemp=initial(1)+r*cos(gtheta);
                        ytemp=initial(2)+r*sin(gtheta);
                        nextelem=obj.Mygeo.findelems([xtemp,ytemp],'in_edge');
                        ratio=ratio+0.01;
                    end
                    [goodelem,pnts,localpnts,Phi]=obj.elementincrement(gtheta,obj.Elemdict(nextelem));
                    % Check if the found element is good
                    if ~goodelem
                        if all(Phi>0)||all(Phi<0)
                            % case a: the crack turns a lot from the current
                            % angle and the above elementincrement method fails
                            % to locate the right nextelem.
                            warning('the nextelem is not correctly determined using the new gtheta');
                            % Approximately use obj.Omega to get the temporary
                            % extention.
                            [~,pnts,localpnts,Phi]=obj.elementincrement(obj.Omega,obj.Elemdict(nextelem));
                            return;
                        else
                            % case b:the found element has little interaction
                            % with the crack--goodelement false
                            initial=pnts(obj.Itip,:);
                            currentelem=obj.Elemdict(nextelem);
                        end
                    end
                end
            else
                nextelem=[];
            end
            obj.NextElem=obj.Elemdict(nextelem);
        end
        
        function [goodelem,pnts,localpnts,Phi]=elementincrement(obj,gtheta,varargin)
            % which is basically a function of OpenGeo fulfilled in Tip
            % 1. Initiate a vector from the current tip
            if ~isempty(varargin)
                elem=varargin{1};
            else
                elem=obj.NextElem;
            end
            %            theta=obj.Growcheck.Growdirection;
            %            gtheta=theta+obj.Omega;
            % Guarantee that the temp segment can cross the whole element
            % and the normal projections of the nodes can lie inside the
            % segment.
            r=3*max(obj.INTELEM.Length);
            xntemp=obj.Xct+r*cos(gtheta);
            yntemp=obj.Yct+r*sin(gtheta);
            % 2. Calculate the signed distance function of nodes in the new
            % element to the new segment, note at least two elements need to
            % update their enrichment after this. (The old tip element and
            % the NextElem)
            tempsegment=[obj.Mygeo.Ends{obj.Itip};xntemp,yntemp];
            if obj.Itip==2
                tempsegment=flipud(tempsegment);
            end
            nodes=elem.NodList;
            plist=[elem.X,elem.Y];
            [xy,distance,~,minind] = distance2curve(tempsegment,plist,'linear');
            minind(minind==3)=2;
            x1=tempsegment(minind,1);
            y1=tempsegment(minind,2);
            x2=tempsegment(minind+1,1);
            y2=tempsegment(minind+1,2);
            xs=xy(:,1);
            ys=xy(:,2);
            x=plist(:,1);
            y=plist(:,2);
            % sign(normal_vec*vec(x-x*))
            signd=sign((y1-y2).*(x-xs)+(x2-x1).*(y-ys));
            Phi=distance.*signd;
            % 3. Return false goodelem if signd has the same sign without
            % runing the intersection method or do it inside the
            % intersection method.
            if all(signd==1)||all(signd==-1)
                goodelem=false;
                pnts=[];
                localpnts=[];
            else
                % 4. Calculate the intersection using the signed distance
                % function, which would be the end point
                [goodelem,pnts,localpnts] = intersection(obj.Mygeo,elem,nodes,Phi);
            end
        end
        
        %% Other function declarations
        initiate(obj);
        [xyl,r]=callocal(obj,x,y);
        % calculate the local coordinates (xl,yl) of a point (x,y), and the
        % distance from the tip
        setradius(obj,varargin);
        calstress_nonlocal(obj);
        realgrow(obj);
    end
end