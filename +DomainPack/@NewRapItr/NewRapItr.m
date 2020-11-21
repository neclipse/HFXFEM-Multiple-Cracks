classdef NewRapItr < handle
   properties
       LinSysCrt        % Class member: Linear system creator and solver
       IncType          % Increment type: 1-load increments; 2-displacement
       IInc=1;          % Index of current load increment
       IItr=0;          % current iteration number
       Steps            % Step sizes for every interval,e.g. [0.1,0.01;1,0.1]
       Tottime          % Total time
       Increments       % Load increments, accumulate percentage of whole load
       MaxInc=1000;      % Control: MAXIMUM number of Increments
       MaxItr=9;        % Control: MAXIMUM Iteration number
       Tolerance=1E-9;  % Tolerance for load vector convergence criteria
       ConvFlag=0;      % Convergence flag
       DivFlag=0;       % Divergence flag
       CutFlag=1;       % Cut flag
       Ratio            % Ratio: the  norm of the residual to the external load vector
       ResMax           % Maximum nodal residual load
       Pincallowed;   % Maximum pressure change allowed per increment, can be seen as an upper limit of pinc
       Pinclimit;  % Threshold to increase the size of increment, can be seen as a lower limit of pinc
       Pincmax;      % Maximu nodal pore pressure increments
       Pratemax;
       
       IntLoadVec 		% Global internal load vector for newton-raphson iteration
       ExtLoadVec		% Global external load vector for newton-raphson iteration
       FullExtLoad      % Full global external load vector;
       IntLoadVecO      % Last Converged internal load vector
       Dt               % Current Time increment
       DtO=0;           % Last time increment
       Newmark          % struct arrary to store parameters for newmark scheme {1*8}
   end 
   properties(Dependent)
       NoInc            % Number of increments
       Timeinc          % Time increments are dependent on the Increments
       Dim              % Dimension of load vector
       ResLoadVec       % Global residual load vector
   end

   methods
       function obj=NewRapItr(inctype,step,tottime,linsyscreater,varargin)
           obj.IncType=inctype;
           obj.Steps=step;
           obj.Tottime=tottime;
           obj.LinSysCrt=linsyscreater;
           obj.IntLoadVec=zeros(obj.Dim,1);
           obj.ExtLoadVec=zeros(obj.Dim,1);   
           obj.IntLoadVecO=zeros(obj.Dim,1);
           if nargin>4 
               obj.MaxInc=varargin{1};
               obj.MaxItr=varargin{2};
               obj.Tolerance=varargin{3};
           end
           if nargin>7
               obj.Pincallowed=varargin{4};
               obj.Pinclimit=varargin{5};
           end
           % Generate the increments array
           obj.incarray;
       end
       function incarray(obj)
           % Generate the increments array
           Incs=[];
           for iint=1:size(obj.Steps,1)
               interval=obj.Steps(iint,1);
               step=obj.Steps(iint,2);
               if iint==1
                   incs=step:step:interval;
               else 
                   incs=obj.Steps(iint-1):step:interval;
               end
               Incs=[Incs,incs];
           end
           obj.Increments=Incs;
           obj.validateinc;
       end
       function validateinc(obj)
           increments=obj.Increments;
           % Guarantee that the last component of the increments array is 1
           if increments(end)<1-1e-4
               increments=[increments,1];
           elseif increments(end)>1
               increments=increments(increments<1);
               increments=[increments,1];
           end
%            if length(increments)>obj.MaxInc
%                warning('Too many increments are involved.');
%            end
           temp=uniquetol(increments);
           temp(1-temp<1e-4)=[];
           obj.Increments=[temp,1];
       end
       
       function NoInc=get.NoInc(obj)
          NoInc=length(obj.Increments); 
       end
       function Dim=get.Dim(obj)
           Dim=obj.LinSysCrt.Drow;
       end
       function Timeinc=get.Timeinc(obj)
          Timeinc=obj.Tottime*obj.Increments; 
       end
       function ResLoadVec=get.ResLoadVec(obj)
          ResLoadVec=obj.ExtLoadVec-obj.IntLoadVec; 
       end
       function set.LinSysCrt(obj,linsyscreater)
           if(isa(linsyscreater,'DomainPack.LinSysCreater'))  
               obj.LinSysCrt=linsyscreater;
           else
               error('Input must be an instance of LinSyscreater');
           end          
       end
       % update the parameters for newmark scheme
       updatenewmark(obj,varargin);
       initialize_v1_E(obj);    % Assume the material at the initial state is pure elastic, not used for elastic analysis
       returning(obj);          % Main plastic returning method: loop over equilibrium equation, not used for elastic XFEM 
       iterating(obj,i,stdpdofs,varargin); % used for elastic
       stagechangeflag=intforcer(obj,varargin);          % Calulate the internal force vector for all elements
       converger(obj,varargin);          % Determine the convergenence flag for the domain
       autoincrem(obj,mode,varargin);  % Automatically adjust the size of the increments
       % IMPORTANT :REMEMBER TO ASSIGN CORRECT INPUT AND OUTPUT NUMBER IN
       % DECLARING THE FUNCTION NAME ONLY
       switching(obj,mode);     % switch between current value and last converged value
       checksize_porepressure( obj,stdpdofs,varargin ); % applying a criterion for automatic step 
       epnotifier(obj);   % not used for pure elastic analysis
   end 

end