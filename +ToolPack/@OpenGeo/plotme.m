function plotme( obj, varargin)
%PLOTME a method of OpenGeo, plot the crack against the mesh
% usage:obj1.plotme(deformflag,crackflag, nodeflag)
% or plotme(obj1,deformflag,crackflag, nodeflag)

% Input Arguments:
% if deformflag == true, then plot the deformed mesh with the displacement
% information, if not, then plot the initial mesh
% if crackflag == true, then plot the crack line against the mesh
% if nodeflag == true, then plot the enriched nodes on the mesh
% plot the mesh, Mesh class has its own method to treat deformflag.
% Default input arguments:
if isempty(varargin)
   deformflag = false;
   crackflag = true;
   nodeflag = true;
   phiflag = false;
else
   deformflag = varargin{1};
   crackflag = varargin{2};
   nodeflag = varargin{3};
   phiflag = varargin{4};
   ux=varargin{5};
   uy=varargin{6};
   scale=varargin{7};
end
if deformflag % plot the crack on the deformed mesh
    VXdeformed=obj.Mesh.VX+ux*scale;
    VYdeformed=obj.Mesh.VY+uy*scale;
    % check if need plot the crack line
%     if crackflag
%         plot(obj.Segments(:,2),obj.Segments(:,3),'-','LineWidth',2);
%         hold on
%     end
    % check if need plot the enriched nodes
    if nodeflag
        bodyx=VXdeformed(obj.Bodynodes);
        bodyy=VYdeformed(obj.Bodynodes);
        tipx=VXdeformed(obj.Rtipnodes);
        tipy=VYdeformed(obj.Rtipnodes);
        plot(bodyx,bodyy,'o');
        plot(tipx,tipy,'s');
    end
    if phiflag
        plus=obj.Phi(:,2)>0;
        allnodes=obj.Phi(:,1);
        plusnode=allnodes(plus);
        xplus=VXdeformed(plusnode);
        yplus=VYdeformed(plusnode);
        minus=obj.Phi(:,2)<0;
        minusnode=allnodes(minus);
        xminus=VXdeformed(minusnode);
        yminus=VYdeformed(minusnode);
        plot(xplus,yplus,'^');
        plot(xminus,yminus,'v');
    end
else            % plot the crack on the undeformed mesh
    if crackflag
        plot(obj.Segments(:,2),obj.Segments(:,3),'-','LineWidth',3);
        hold on
    end
    % check if need plot the enriched nodes
    if nodeflag
        bodyx=obj.Mesh.VX(obj.Bodynodes);
        bodyy=obj.Mesh.VY(obj.Bodynodes);
        tipx=obj.Mesh.VX(obj.Rtipnodes);
        tipy=obj.Mesh.VY(obj.Rtipnodes);
        plot(bodyx,bodyy,'o');
        plot(tipx,tipy,'s');
    end
    if phiflag
        plus=obj.Phi(:,2)>0;
        allnodes=obj.Phi(:,1);
        plusnode=allnodes(plus);
        xplus=obj.Mesh.VX(plusnode);
        yplus=obj.Mesh.VY(plusnode);
        minus=obj.Phi(:,2)<0;
        minusnode=allnodes(minus);
        xminus=obj.Mesh.VX(minusnode);
        yminus=obj.Mesh.VY(minusnode);
        plot(xplus,yplus,'^');
        plot(xminus,yminus,'v');
    end
end
end

