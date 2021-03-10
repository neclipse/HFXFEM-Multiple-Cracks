function autoincrem( obj, mode, varargin )
%AUTOINCREM Method of class NewRapItr, accumulated increments
% Modified on 07/12/2019 to accomodate unstable crack growth by
% inserting a series of tiny time steps after the current increment
% and reduce the following time step by dratio.
%% Input parser
p=inputParser;
addOptional(p,'dratio',0.25);
addOptional(p,'iratio',1.5);
addOptional(p,'inclist',[]);
addOptional(p,'minimalinc',0.0001);
addOptional(p,'allwoedsteps',3);
parse(p,varargin{:});
dratio=p.Results.dratio;
iratio=p.Results.iratio;
inclist=p.Results.inclist;
minimalinc=p.Results.minimalinc;
allowedsteps=p.Results.allwoedsteps;
iinc=obj.IInc;
switch mode
    %  Auto Cut the current load increment factor by 50 percent
    case 1
        fprintf('Increment size is cut at No.%d increment\n',obj.IInc);
        if iinc==1
            newincrem=dratio*obj.Increments(iinc);
            %             obj.Increments=newincrem:newincrem:1;
            increments=[newincrem,obj.Increments(2:end)];
        elseif iinc>1
%             %% Approach 1: May introduce too many increments
%             step=dratio*(obj.Increments(iinc)-obj.Increments(iinc-1));
%             newincrem=obj.Increments(iinc-1)+step;
%             front=obj.Increments(1:iinc-1);
%             tail=obj.Increments(iinc:end);
%             increments=[front,newincrem,tail];
%             % To find the checkpoint within obj.Steps
%             temp1=obj.Steps(:,1);
%             % BUG: SHOULD NOT USE TEMP1(TEMP1>0) TO FIND THE MINIMUM POSITIVE 08092019
%             checkpoint=min(temp1(temp1>=newincrem));
%             nextstage=newincrem:step:checkpoint;
%             tail= obj.Increments(obj.Increments>checkpoint);          
%             increments=[front,nextstage,tail];
            %% Approach 2: Did not significantly affect the number steps afterawards
            front=obj.Increments(1:iinc-1);
            current=obj.Increments(iinc-1:end-1);
            back=obj.Increments(iinc:end);
            steps=dratio*(back-current);
            newincrems=current+steps;
            increments=[front,newincrems];
        end
        % Auto elongate the subsequent load incremnt factor by 25 percent
    case 2
%         fprintf('Increment size is elongated after No.%d increment\n',obj.IInc);
        step=iratio*(obj.Increments(iinc+1)-obj.Increments(iinc));
        newincrem=obj.Increments(iinc)+step;
        front=obj.Increments(1:iinc);
        %         back=obj.Increments(iinc+2:end);
        if newincrem<1
            if iinc+2<obj.NoInc
                if newincrem+step<obj.Increments(iinc+2)
                    back=obj.Increments(iinc+2:end);
                    obj.Increments=[front,newincrem,back]; % to take account for the change in obj.Steps(:,1)
                    return;
                end
            end
            istep=find(newincrem<obj.Steps(:,1),1);
            back=newincrem:step:obj.Steps(istep,1);
            [~,iincrement]=min(abs(obj.Increments-obj.Steps(istep,1)));
            rest=obj.Increments(iincrement+1:end);
            increments=[front,back,rest];
        end
    case 3 % for unstablegrow insert several small time steps
        fprintf('Increment size is cut after No.%d increment to allow unstable crack growth\n',obj.IInc);
        front=obj.Increments(1:iinc-1);
        insertion=obj.Increments(iinc):minimalinc:(obj.Increments(iinc)+minimalinc*allowedsteps);
        current=obj.Increments(iinc:end-1);
        back=obj.Increments(iinc+1:end);
        steps=dratio*(back-current);
        newincrems=current+steps;
        % Need to delete the first newincrems as it may be smaller than the
        % current increment so that the current increment will be repeated
        % after sorting. 07142019
        %         newincrems(1)=[];
        increments=[front,insertion,newincrems];
end
increments=[increments,inclist];
obj.Increments=increments;
obj.validateinc;
end


