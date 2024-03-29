function obj=preparing(obj,X,Y, EnrichNum)
% method of GaussPnt_Cohesive class: TO WATCH. DOES THE PREPARING METHOD
% NEED TO BE CHANGED FOR THE LINE GAUSSIAN POINTS. NOTE THIS METHOD WAS
% ORIGINALLY DEVELOPED FOR THE 2D DOMAIN GAUSSIAN POINTS. 10/20/20

% It provides Jacobian and the derivative of shape functions over global coordinates.
% It also assemble the B matrix for the gaussian point.
% It also gives the coordinates of the gausssian point in global
% coordinate system
    nnodes=obj.NNodes;
    xi=obj.Xi;
    eta=obj.Eta;
    X=X(1:nnodes);
    Y=Y(1:nnodes);
    N_x=zeros(1,nnodes);
    N_y=zeros(1,nnodes);
    [N,N_xi,N_eta]=shape(nnodes,xi,eta);
    x_xi=N_xi*X;
    x_eta=N_eta*X;
    y_xi=N_xi*Y;
    y_eta=N_eta*Y; 
    J=[x_xi y_xi;x_eta y_eta];
    dJ=det(J);% determinant of J                                                                                           
    for k=1:nnodes
        temp=J\[N_xi(k);N_eta(k)];
        N_x(k)=temp(1);
        N_y(k)=temp(2);
    end
    % Create Nu matrix
    Nu=zeros(2,2*nnodes);
    % PREALLOCATE NUENR MATRIX 10/19/20
    obj.Nuenr=zeros(2,2*EnrichNum*nnodes);
    obj.Nuenrplus=obj.Nuenr;
    obj.Nuenrminus=obj.Nuenr;
    Nu(1,1:2:2*nnodes-1)=N;
    Nu(2,2:2:2*nnodes)=N;
    % Create B matrix
    B=zeros(4,2*nnodes); % plane strain problem,dimension(S)=[4,2]
    % PREALLOCATE BMATENR
    obj.Bmatenr=zeros(4,2*EnrichNum*nnodes);
    % Assemble B matrix row by row
    B(1,1:2:2*nnodes-1)=N_x;
    B(2,2:2:2*nnodes)=N_y;
    B(3,1:2:2*nnodes-1)=N_y;
    B(3,2:2:2*nnodes)=N_x;
    % Create DNp matrix
    obj.Np=N;
    obj.Npenr=zeros(size(N,1),size(N,2)*EnrichNum);
    obj.Nu=Nu;
    obj.DNp=[N_x;N_y];
    obj.DNpenr=zeros(size(obj.DNp,1),size(obj.DNp,2)*EnrichNum);
    obj.Bmat=B;
    obj.DetJ=dJ;
    % locate gaussian point in global cartesian coordinate system
    obj.X=N*X;              
    obj.Y=N*Y;
end
