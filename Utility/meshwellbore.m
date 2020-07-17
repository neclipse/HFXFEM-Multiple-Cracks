function [ node,element,sznode,szelement ] = meshwellbore( re,rw,partition,tanstep,radstep1,radstep2,bias1,bias2 )
%MESH mesh preparation for nodes info and element info
%% Partition the mesh domain by the partition ratio
length=re-rw;
length1=(1-partition)*length;                       % The radial length of the part1, inner part
length2=partition*length;                           % The raidal length of the part2
ri1=rw;                                             % The radial coordinate of the starting point of part2 
ri2=ri1+length1;                                    % The radial coordinate of the starting point of part2 
%% Discretize each partition
omega=pi/2/tanstep;                                 % The tangential angle step size
tannode=tanstep+1;                                  % The number of nodes along tangential direction
sznode=(radstep1+radstep2+1)*tannode;                        % The quantity of nodes
szelement=(radstep1+radstep2)*tanstep;                         % The quantity of elements
node=zeros(sznode,2);                               % prepare a zero matrix to store info of nodes including sequence No.,x and y
element=zeros(szelement,4);                         % prepare a zero matrix to store info of elements
%%--To assign value to nodes
%-- Prat1
% calculate the first element length
q1=bias1^(1/(radstep1-1));                          % The scale factor of part1
li1=length1*(1-q1)/(1-q1^radstep1);                   % The length of the smallest element
k=1;
r=ri1;
for i=1:radstep1+1
    if i>1
    lnew=li1*(q1^(i-2));
    r=r+lnew;
    end
    for j=1:tannode
        node(k,1)=r*cos((j-1)*omega);% x coordinate value
        node(k,2)=r*sin((j-1)*omega);% y coordinate value
        k=k+1;
    end
end
%-- part2
% calculate the first element length
q2=bias2^(1/(radstep2-1));                          % The scale factor of part1
li2=length2*(1-q2)/(1-q2^radstep2);                   % The length of the smallest element
r=ri2;
for i=1:radstep2
    lnew=li2*(q2^(i-1));
    r=r+lnew;
    for j=1:tannode
        node(k,1)=r*cos((j-1)*omega);% x coordinate value
        node(k,2)=r*sin((j-1)*omega);% y coordinate value
        k=k+1;
    end
end

%%--To assign value to elements
k=1;
for i=1:(radstep1+radstep2)
    for j=1:tanstep
        initial=(i-1)*tannode+j;
        %one line stores the sequence No. of all four nodes in one element in anticlockwise order
        element(k,:)=[initial,initial+tannode,initial+tannode+1,initial+1]; 
        k=k+1;
    end
end
end

