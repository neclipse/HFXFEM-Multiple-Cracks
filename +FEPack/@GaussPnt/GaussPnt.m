classdef GaussPnt < handle
    properties
    Xi          % s-coordinate in local parent coordinate system
	Eta         % t-coordinate in local parent coordinate system
	H           % integration weight
    X           % x-coordinates in global cartesian coordinate system
    Y           % y-coordinates in global cartesian coordinate system
    NNodes      % number of Nodes involved in the calculation of gaussian points
    Bmat        % B matrix for this gaussian point
    DetJ        % Determinant of the Jacobian matrix
    GBINP       % Global inputs that required for matsu
    end
    
    methods
		function obj=GaussPnt(xi,eta,h,nnodes,gbinp)
            if nargin>0
			obj.Xi=xi;
			obj.Eta=eta;
			obj.H=h;
            obj.NNodes=nnodes;             % number of nodes involved in the calculation of gaussian points
            obj.GBINP=gbinp;
            end
        end
        function initial(obj,inistress,Delastic)% Initialize the elastic strain of the prestressed material
            obj.ElaStnO=Delastic\inistress;
        end
        preparing(obj,X,Y); % Calculate X-,Y- coordinates, B matrix, dJ, and other derivative for the left handside of the linear system
        switching(obj,mode);
    end
    methods(Abstract)
		matsu(obj);
		matct(obj);
	end
end