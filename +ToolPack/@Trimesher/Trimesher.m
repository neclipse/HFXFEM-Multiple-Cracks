classdef Trimesher < ToolPack.Mesher2d
    
    methods
        %Constructor
        function obj=Trimesher(fd,fh,h0,bbox,pfix)
            %        fd        %Distance function d(x,y)
            %        fh        %Scaled edge length function h(x,y)
            %        h0        %Initial edge length
            %        bbox      %Bounding box [xmin,ymin; xmax,ymax]
            %        pfix      %Fixed node positions (NFIXx2)
            if nargin==0
                error('Domain geometry input missing')
            else
            [obj.p, obj.t]=distmesh2d(fd,fh,h0,bbox,pfix);
            end
            
        end
        
        %Validate the quality of the mesh
        function validate(obj)
           obj.q=MeshQuality(obj.EToV,obj.VX,obj.VY);
           maxit=10; 
           while min(obj.q)<0.5
               fprintf('Bad elements dectected,smooth function will be called');
               obj.smooth(maxit,0.05)
               maxit=maxit*log(maxit);
               obj.q=MeshQuality(obj.EToV,obj.VX,obj.VY);
           end
           fprintf('The mesh of good quality with worst element quality of %.4f\n',min(obj.q));
        end
        
        % Smooth the mesh if quality is not satisfied(<0.5)
        function smooth(obj,maxit,tol)
            [obj.p, obj.t]=smoothmesh(obj.p,obj.t,maxit,tol);
            disp('meshsmooth function called')
            % maxit : maximum number of iterations
            % tol:Convergence tolerance (Percentage change in edge length
            % must be less than TOL).
        end
        
        % Reorder the mesh nodes according to the convention
        function reorder(obj)
            obj.t=Reorder(obj.t,obj.VX,obj.VY);
            disp('Reorder called')
        end

       function material(obj)
           obj.MT=zeros(size(obj.p,1),3);
           for i=1:obj.Totelem
               xc=sum(obj.VX(obj.EToV(i,:)))/4; % divided by 4 is for the four noded element
               yc=sum(obj.VY(obj.EToV(i,:)))/4;
               obj.MT(i,:)=[1,xc,yc]; % '1' is the material property type, which is set for all elements temporarily
           end

       end
       
       function nodetype(obj)
          obj.VT=ones(size(obj.p,1),1);
       end
    end
end
