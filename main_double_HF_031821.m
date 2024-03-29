% This is the main script of the proposed enriched approach for the
% stationary crack. 
%The example would be a infinitely large plate with a center crack, due to 
% the symmetry about the center axis, only half mesh is used. 
%%  Importing packages
% -- Initialize the work space
clear; clc;
% add the home folder to the matlab search path. Already saved for future
% session, so I commented this line.
% addpath('C:\Users\chuan\Documents\HFXFEM-Verified-Github\HFXFEM-Verified-Singlecrack');
% Import class packages 
import ToolPack.*;                          % Complimentary tool classes, like Preprocessor
import FEPack.*;                            % Elementary level classes, like element class, node class and gausspnt class
import DomainPack.*;                        % Domain level classes, like Domain, LinSysCrt, NewRapItr
import EnrichPack.*;
% addpath(genpath('.\'));
% Set up the model
%% --Preparation for data
%% --Geometry, loads, BCs settings in Preprocessor
% The plate has a width of lw and a height of lh, and the center linear
% crack has length of 40cm.
lw=45;                                      % The width of the plate
lh=60;                                      % The height of the plate
lc=0.08;                                    % The center crack length 
ls=8;                                       % The spacing between two HF
% four boundary handles (note that these handles should be adjusted based
% on where the origin is located on the mesh.)
fb1= @(p) p(:,2);                    % The bottom side
fb2= @(p) lw-p(:,1);                    % The right side
fb3= @(p) lh-p(:,2);                    % The top side
fb4= @(p) p(:,1);                    % The left side (symmetry axis)
fb={fb1,fb2,fb3,fb4};
% injection
qline=-0.011248593925759;              % Unit flux rate, which is volumetric rate/area: [m/s], negative for source
% the sectional area of the perforated wellbore: wellbore 3.75inch diameter (0.085725 m) ,
% 0.35 inch (0.00889m ) peforation hole; gives:
q=-0.0005;            % flow rate: [m^2/s] -1e-3
% It is kind of realistic as it corresponds to 0.079 bpm per
% perforation hole. In real case, there are more than 5 perforation
% holes per cluster, and more than 5 clusters per stage. The total
% injection rate for a stage ranges from 10 bpm to 60 bpm. Also, it can
% be increased to 1e-3.
%     q=-0.0005;
% meshing
% meshnode=readmatrix('nodefile_01122020_carrier.txt','Delimiter',',');
% meshelement=readmatrix('elemfile_01122020_carrier.txt','Delimiter',',');
meshnode=readmatrix('nodefile_Double_HF_0310.txt','Delimiter',',');
meshelement=readmatrix('elemfile_Double_HF_0310.txt','Delimiter',',');
meshnode=meshnode(:,2:3);
meshelement=meshelement(:,2:5);
plate=Quadmesher(meshnode,meshelement);
% plate.plotmesh;
%% Inputs and Parfor setting
% % To display the progress
% proT=uitable('ColumnName','numbered');
% proT.Units='Normalized';
% proT.Position=[.1,.1,.8,.8];
workers=int8(4);
cases=1;
caseids=1;
% proT.Data(1,1:cases)=caesids;
% parprofile=Par(cases); % profile the time for each head
% parfor(icase=1:cases,workers) % request 4 workers
for icase=1:cases
%     Par.tic;
    fprintf('cases_%d is running\n',caseids(icase));
    GBINP=assembleglobalinputs(caseids(icase));
%% --  Step 1
    %% set up the boundary condtion

    % ----Special nodes
    BCTable_node=[];        % For Dirichelet boundary condition on sparse nodes
    % Do not use QTable_node unless this is injection to non-crack area.
    % 09082019
    QTable_node=[];
%     QTable_node=[[1,4],{qpoints},q];   % For line source/link in arbitrary location
    % ---- Comprehensive boundary conditions for the external boundary
    % describe by fb. fixed bottom boundary and permeable, tensile traction
    % applied to the top and injection is also applied through the top. The
    % rest two boundaries are impervious and free from deformation in y
    % direction.
    BCmat_line=[1,1,1,0;1,1,1,0;1,1,1,0;1,0,0,0];% The description is in the Precessor class.
%     BCmat_line=[1,1,1,0;1,1,1,0;1,0,1,4;1,0,0,0];% The description is in the Precessor class.
    [ir,ic]=find(BCmat_line==1);
    psdboundind=[ir,ic];    % Dirichlet boundary, in this case all bounaries are fixed
    imbalancedbound=[];     % Neumann boundary is the 3rd boundary, top edge, must specify if existing  
    % THE SPECIFIED BOUNDARIES WILL BE USED IN LINSYS.INITIALRHS AND THE
    % INSITU STRESS RESULTANT LOAD VECTOR WILL BE CALCULATED. THE CURLOAD
    % SHOULD BE SPECIFIED AND THE BCMAT_LINE SHOULD CHANGE FROM '1' TO '4'.
    predof=[0,0,0;0,0,0;0,0,0;0,0,0]; % Predefined dofs corresponding to the '1's in the BCmat_line;
%     curloads1m=[0,3.7e-3];   % [Tx,Ty] are the current traction to the top side, GPa 
%     curloads3m=[0,-5.2e-3];   % [Tx,Ty] are the current traction to the top
%     side, GPa. 
    curloads3m=[];
%     curloads3q=q;       % OUTWARD POSITIVE, INJECTION NEGATIVE.
    curloads={[],[];[],[];curloads3m,[];[],[]};
    %% Preprocessing
    preprocessor1=ToolPack.Preprocessor(4,fb,BCmat_line,BCTable_node,predof,curloads,psdboundind,imbalancedbound,QTable_node);
    preprocessor1.Mesher=plate;
    preprocessor1.preprocess;
    % Set the domain with UP element, poroelastic material, precribed displacements to boundaries
    Step1=DomainPack.Domain(4,4);
    Step1.GBINP=GBINP;
    Step1.Preprocess=preprocessor1;           % Assign the defined preprocessor to a property of the Wb object
    % Meshing: The generated coordinates are defined as the initial location of
    % every nodes, further deformation is defined with respect to this set of coordinates.
    Step1=Step1.assignmesh;
    Step1=Step1.updatedofarray;
    Step1=Step1.crtlinsys;                   % create the linsystem here because enrichitem may be void (03032019,06232019)
    %% Setting up the crack profile, consider to add multiple cracks
    mesh=Step1.Preprocess.Mesher;
    elemdict=Step1.ElemDict;
    nodedict=Step1.NodeDict;
    bdls=Step1.Preprocess.Disthandle;
    % set crack geometry using function handles
%     crackfunhandle=@(x) 0.5*x;
%     xlimits=[-0.10,0.10];
%     des={crackfunhandle,1,xlimits};
%     crack1=ToolPack.OpenGeo(1,mesh,bdls,nodedict,elemdict,2,des,10);
   % set crack geometry using segment points for Two HFs
    segments1=[1,0,33.9996;2,lc,33.9996];      % The mesh file is a little bit shifted, to adjust.         
    segments2=[1,0,26.0001;2,lc,26.0001];  % crack segments [n,x,y]
    % x should always start from left to right. 
    % Enforced this condition in opengeo.discretize 01/02/21
    % set crack geometry for two NFs (3m long=1+2)
    theta1=30*pi/180;
    theta2=15*pi/180;
    theta3=60*pi/180;
    segments3=[1,1.5-0.1*cos(theta3),lh/2+ls/2-0.1*sin(theta3);2,1.5+0.2*cos(theta3),lh/2+ls/2+0.2*sin(theta3)]; 
    segments4=[1,0.6-0.05*cos(theta2),lh/2-ls/2-0.05*sin(theta2);2,0.6+0.2*cos(theta2),lh/2-ls/2+0.2*sin(theta2)]; 
    segments5=[1,2.5-0.5*cos(theta1),lh/2+ls/2-0.5*sin(theta1);2,2.5+1*cos(theta1),lh/2+ls/2+1*sin(theta1)];
    HF1=ToolPack.OpenGeo(1,mesh,bdls,nodedict,elemdict,1,segments1,10); % The HF1
    HF2=ToolPack.OpenGeo(2,mesh,bdls,nodedict,elemdict,1,segments2,10); % The HF2
    NF1=ToolPack.OpenGeo(3,mesh,bdls,nodedict,elemdict,1,segments3,10); % The NF1
    NF2=ToolPack.OpenGeo(4,mesh,bdls,nodedict,elemdict,1,segments4,10); % The NF2
    NF3=ToolPack.OpenGeo(5,mesh,bdls,nodedict,elemdict,1,segments5,10); % The NF3
    crackdict=[HF1,HF2,NF1,NF2,NF3];
    injectionpoint1=[0,lh/2+ls/2]; % needed for opengeo.findblending.
    injectionpoint2=[0,lh/2-ls/2]; % needed for opengeo.findblending.
    % visual check of the cracks and the nodes detection.
%     mesh.plotmesh;
%     hold on;
    for icrack=1:length(crackdict)
        crackdict(icrack).initiate;
%         crackdict(icrack).plotme;
    end
    % set EnrCrackBody using the initial crack info
    % Initialmode 1: perforated, completely traction free; 
    % 2: "smeared crack" or cemented crack, seemingly continuum
    % 3:existing fracture, start with compressive mode
    % 4: existing fracture, start with tensile mode.
    % 5: newly propagated fracture started from tensile mode.
    InitialMode1=1; % 1:perforated for HF
    InitialMode2=1; 
    InitialMode3=1; % 2: smeared for NF. 
    InitialMode4=1;
    InitialMode5=1;
    cohesivetype='unified';
    % This alpha is used to initiate the initial traction and crack opening
    % for existing open crack with cohesive traction, implemented in 
    % GaussPnt_Cohesive class. In fact, this is no longer useful as the
    % existing open crack is usually modeled as perforated 'true' and there
    % is no remaining cohesion to initiate with. For the initially existing
    % weak crack (modeled as smeared rock matrix), the angle is not useful
    % either. Because the traction will be initiated from the stress at the
    % linegauss points along the predefined crack path. The linegauss and
    % subdomain gaussian points could be prepared even when they the element
    % are not enriched yet. The enrichment will be hold until the
    % calculated cohesive traction exceeds the threshold. 11/04/2020
    % This approach although more complicated, would be more robust than
    % the current handling of contact modes.
%     Alpha1=pi/2;     %pi/2-atan(0.5) the angle between the tensile force and crack plane
%     Alpha2=pi/2; % not necessary as perforated is true. 
    encrack1=EnrichPack.EnrCrackBody('crackbody',elemdict,nodedict,HF1,InitialMode1,cohesivetype);
    encrack2=EnrichPack.EnrCrackBody('crackbody',elemdict,nodedict,HF2,InitialMode2,cohesivetype);
    encrack3=EnrichPack.EnrCrackBody('crackbody',elemdict,nodedict,NF1,InitialMode3,cohesivetype);
    encrack4=EnrichPack.EnrCrackBody('crackbody',elemdict,nodedict,NF2,InitialMode4,cohesivetype);
    encrack5=EnrichPack.EnrCrackBody('crackbody',elemdict,nodedict,NF3,InitialMode5,cohesivetype);
    encrack1.Qtable=[encrack1.Id,q]; % for edge crack, it is okay to ignore the injection point.
    encrack2.Qtable=[encrack2.Id,q]; % for edge crack, it is okay to ignore the injection point.
    Step1.EnrichItems=[encrack1,encrack2,encrack3,encrack4,encrack5];           
    %% Start the Newton-Raphson iterative analysis
    % ---- Newton-Raphson Iterator
%     step=[0.005,0.0001;0.2,0.006;1,0.01];          % dimensionless increment size
%     step=[0.005,0.0002;0.3,0.003;1,0.009];          % dimensionless increment size
%     step=[0.01,0.00015;0.1,0.0008;0.4,0.002;0.8,0.004;1,0.008];  % dimensionless increment size
    step=[0.01,0.0002;0.1,0.001;0.4,0.004;0.8,0.008;1,0.01];  % dimensionless increment size
    tottime=10;                                 % total time
    inctype=1;                                  % inctype: 1-load increments; 2- displacement increments
    % The following three parameters are optional setting to control the speed
    % and the accuracy3 of the Newton-Raphson algorithm. If one is not sure the
    % impact of the setting, one can leave them blank.
    maxinc=2000;
    maxitr=15; % set to an odd number upon issue # 35.
    tol=1E-7;
    pincallowed=[];      % upper limit for pinc, may cause continuous increment cut if too small
    pinclimit=0.00001;       % a threshold to increase increment size
    newraph1=DomainPack.NewRapItr(inctype,step,tottime,Step1.LinSysCrt,maxinc,maxitr,tol,pincallowed,pinclimit);
%     newraph1=NewRapItr(inctype,step,tottime,Step1.LinSysCrt,maxinc,maxitr,tol);
    Step1.NewtonRaphson=newraph1;
%     % To store the results in postprocessor
    postdict(1,maxinc)=ToolPack.Postprocessor();
    Step1.Postprocess=postdict;
    % Running the simulation
    % -- Begin loop over load increments
%     tic;                                        % set a timer
    % Wb.NewtonRaphson.returning;
    % savemode defines the how often the results are saved: 1--means every increment; 2--means every two; n--every n increments
%     savemode=1;
%     Step1=Step1.running(postdict,savemode);
    savemode=2;
    saveinc=int8(10);
    Step1=Step1.running(postdict,savemode,'interval',saveinc);
%     savemode=3;
%     steplist1=[0.01,0.1,1];
%     Step1=Step1.running(postdict,savemode,'inclist',steplist1);
%     save('propagating_comparison_1105.mat','Step1');
%     toc;
    postdict=Step1.Postprocess;
    filename=strcat('C:\Users\chuan\Google Drive\Exciting Research\Writings\Efficient HM-XFEM model with complex fracture network\Results\Enhanced rerun-0321\case_',num2str(caseids(icase)),'_newangles5.mat');
    m=matfile(filename,'writable',true);
    m.postdict=postdict;
    m.GBINP=GBINP;
    fprintf('cases_%d is finished\n',caseids(icase));
end
% poolobj=gcp('nocreate');
% delete(poolobj);
%% Plotting
%example of postprocessing
% export_fig 'Organized results\Stationary_center_crack_0410\Deformed mesh plot with center crack.tif' -m3.125 -transparent
% mesh.plotmesh;
% % 
figure()
hold on;
for icrack=1:length(crackdict)
crackdict(icrack).plotme; % deform, crack, node, phi,ux,uy,scale
end
ax=gca;
fs=16;
titlestr=strcat('The final fracture network for case-3 by 60s injection');
% title(titlestr,'FontSize',fs);
axis('equal')
xlabel('X axis')
ylabel('Y axis')
xlim([0,10])
ylim([24,38])
% Set x and y font sizes.
ax.XAxis.FontSize = fs-2;
ax.YAxis.FontSize = fs-2;
legend('HF1','HF2','NF1','NF2','NF3','FontSize',fs-2)