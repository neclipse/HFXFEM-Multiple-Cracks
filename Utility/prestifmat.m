function [dJ,N,N_x, N_y,B]=prestifmat(X,Y,xi,eta,nnodes)
        % This function prepare for the stiffness matrix.
        % It provides Jacobian and the derivative of shape functions over global coordinates.
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
        % Create B matrix
        B=zeros(4,2*nnodes); 
        % Assemble B matrix row by row
        B(1,1:2:2*nnodes-1)=N_x;
        B(2,2:2:2*nnodes)=N_y;
        B(3,1:2:2*nnodes-1)=N_y;
        B(3,2:2:2*nnodes)=N_x;
end