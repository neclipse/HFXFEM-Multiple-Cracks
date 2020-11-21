% This is the main script of the proposed enriched approach for the
% stationary crack. 
%The example would be a infinitely large plate with a center crack, due to 
% the symmetry about the center axis, only half mesh is used. 
%%  Importing packages
% -- Initialize the work space
clear; clc;
% Import class packages 
import ToolPack.*;                          % Complimentary tool classes, like Preprocessor
import FEPack.*;                            % Elementary level classes, like element class, node class and gausspnt class
import DomainPack.*;                        % Domain level classes, like Domain, LinSysCrt, NewRapItr
import EnrichPack.*;
addpath(genpath('.\Utility'));
% Set up the model
%% --Preparation for data
%% --Geometry, loads, BCs settings in Preprocessor
% The plate has a width of lw and a height of lh, and the center linear
% crack has length of 40cm.
lw=45;                                      % The width of the plate
lh=60;                                      % The height of the plate
lc=0.05;                                    % The center crack length     
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
meshnode=readmatrix('nodefile_Dontsov_medium.txt','Delimiter',',');
meshelement=readmatrix('elemfile_Dontsov_medium.txt','Delimiter',',');
meshnode=meshnode(:,2:3);
meshelement=meshelement(:,2:5);
plate=Quadmesher(meshnode,meshelement);
% plate.plotmesh;
%% Inputs and Parfor setting
% % To display the progress
% proT=uitable('ColumnName','numbered');
% proT.Units='Normalized';
% proT.Position=[.1,.1,.8,.8];
% workers=int8(4);
% cases=12;
% caseids=71:82;
% % proT.Data(1,1:cases)=caesids;
% % parprofile=Par(cases); % profile the time for each head
% parfor(icase=1:cases,workers) % request 4 workers
%     Par.tic;
%     fprintf('cases_%d is running\n',caseids(icase));
    GBINP=assembleglobalinputs(1);
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
%     BCmat_line=[0,1,1,0;0,0,1,0;0,0,1,4;1,0,0,0];% The description is in the Precessor class.
    [ir,ic]=find(BCmat_line==1);
    psdboundind=[ir,ic];    % Dirichlet boundary, in this case all bounaries are fixed
    imbalancedbound=[];     % Neumann boundary is the 3rd boundary, top edge, must specify if existing  
    % THE SPECIFIED BOUNDARIES WILL BE USED IN LINSYS.INTIALRHS AND THE
    % INSITU STRESS RESULTANT LOAD VECTOR WILL BE CALCULATED. THE CURLOAD
    % SHOULD BE SPECIFIED AND THE BCMAT_LINE SHOULD CHANGE FROM '1' TO '4'.
    predof=[0,0,0;0,0,0;0,0,0;0,0,0]; % Predefined dofs corresponding to the '1's in the BCmat_line;
%     curloads1m=[0,3.7e-3];   % [Tx,Ty] are the current traction to the top side, GPa 
%     curloads3m=[0,-2e-3];   % [Tx,Ty] are the current traction to the top
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
   % set crack geometry using segment points
    segments=[1,0,lh/2;2,lc,lh/2];                      % crack segments [n,x,y]
    des=segments;
    crack1=ToolPack.OpenGeo(1,mesh,bdls,nodedict,elemdict,1,des,10);
    injectionpoint=[0,lh/2];
    crack1.initiate_2nd(injectionpoint);
    Perforated=true;
    Cohesive='unified';
    Alpha=pi/2;     %pi/2-atan(0.5) the angle between the tensile force and crack plane
    % This alpha is used to initiate the initial traction and crack opening
    % for existing crack, implemented in GaussPnt_Cohesive class.
    encrack=EnrichPack.EnrCrackBody(1,'crackbody',elemdict,nodedict,crack1,Perforated,Cohesive,Alpha);
    % for injection, both qext_enr and qext_std
    crackqtable=[encrack.Id,q]; % for edge crack, it is okay to ignore the injection point.
    encrack.Qtable=crackqtable;
    enfh=EnrichPack.EnFHeaviside(1,crack1.Minelength,crack1.Phi); % smoothed type 1 heaviside is required
    enfrd=EnrichPack.EnFRidge(crack1.Phi);
    % May add another junction enrichment function for minor crack or major
    % crack or separately or both?
    encrack.Myenfs={enfh,enfrd};
    encrack.initial_enrich;
    % Explicitly assign the propagation check and direction 
%     encrack.Mytips(1).GrowCheck=ToolPack.Maxpscheck(threshold);
%     encrack.Mytips(1).Growdirection=ToolPack.Maxpsdirection(encrack.Mytips(1).Omega);
%     encrack.Mytips(2).GrowCheck=ToolPack.Maxpscheck(threshold);
%     encrack.Mytips(2).Growdirection=ToolPack.Maxpsdirection(encrack.Mytips(2).Omega);
    Step1.EnrichItems={encrack}; % need to be updated when new iterms are added 090920
%     mesh.plotmesh;hold on; crack1.plotme;
    % -- Creating a sparse linear system upon updating the dofarray             
    Step1=Step1.updatedofarray_enriched;                        % update the linsystem inside
    %% Start the Newton-Raphson iterative analysis
    % ---- Newton-Raphson Iterator
%     step=[0.005,0.0001;0.2,0.006;1,0.01];          % dimensionless increment size
%     step=[0.005,0.0002;0.3,0.003;1,0.009];          % dimensionless increment size
    step=[0.004,0.0001;0.2,0.003;0.75,0.008;1,0.01];  % dimensionless increment size
    tottime=15;                                 % total time
    inctype=1;                                  % inctype: 1-load increments; 2- displacement increments
    % The following three parameters are optional setting to control the speed
    % and the accuracy3 of the Newton-Raphson algorithm. If one is not sure the
    % impact of the setting, one can leave them blank.
    maxinc=1000;
    maxitr=12;
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
    

     %% Storing the Numerical solutions, Only needed for parfor 090920
%     startpoint=1;
%     endpoint=length(postdict);
%     % endpoint=61;
%     NT=transpose([postdict(startpoint:endpoint).Inc]);
%     NL=zeros(size(NT));
%     NCMOD=NL;
%     NCMP=NL;
%     CKV=NL;
%     LKV=NL;
%     % samplenodes=[445,4836,1615,1635,1645,1705,1755,1865,1915,1945]; % mesh
%     % Dontsov 0403
%     % samplenodes=[9333,9317,9297,9277,9257,9237,9217,9207]; % after 02/24/20
%     % samplenodes=[9337,9317,9297,9277,9257,9237,9217];   % prior to 02/24/20
%     % LKF=zeros(length(NL),length(samplenodes));
%     % CL=zeros(1,length(samplenodes));
%     for i=startpoint:endpoint
%         mouthind=1;
%         NL(i-startpoint+1)=postdict(i).EnrichItems{1}.Length;
%         NCMOD(i-startpoint+1)=postdict(i).EnrichItems{1}.Aperture(mouthind);
%         NCMP(i-startpoint+1)=postdict(i).EnrichItems{1}.Pfrack(mouthind);
%         LKV(i-startpoint+1)=postdict(i).EnrichItems{1}.LeakoffVolume;
%         CKV(i-startpoint+1)=postdict(i).EnrichItems{1}.CrackVolume;
%         %     LKF(i-startpoint+1,:)=postdict(i).EnrichItems{1}.LeakoffFlux;
%     end
%     NL=(NL-0.05);  % 0.05 is the initial perforated length
%     NCMOD=(NCMOD-GBINP.perfaperture)*1e3; % mm
%     NCMP=NCMP*1e3; %MPa
%     Results=[NT,NL,NCMOD,NCMP,CKV,LKV];
%     filename=strcat('case_',num2str(caseids(icase)),'.mat');
%     m=matfile(filename,'writable',true);
%     m.postdict=postdict;
%     m.GBINP=GBINP;
%     m.Results=Results;
% %     parprofile(icase)=Par.toc;
%     fprintf('cases_%d is finished\n',caseids(icase));
%     m.time=time
% end
% stop(parprofile);
% plot(parprofile);
% poolobj=gcp('nocreate');
% delete(poolobj);
%% Plotting
%example of postprocessing
% export_fig 'Organized results\Stationary_center_crack_0410\Deformed mesh plot with center crack.tif' -m3.125 -transparent