function snapshot( obj, crackid, variable, timesteps, varargin)
%SNAPSHOT show snapshots of results from the postprocess
% obj      : is the domain object
% variable : is the primary variable to be visualized: p, ux, etc
% especially useful for the enrichitem to show the evolution of crack
% timesteps: the required timesteps at which the snapshots are to be
% generated.
p=inputParser;
addOptional(p,'component',1);
addOptional(p,'timeincs',[]);
addOptional(p,'basecoordinate',0);
parse(p,varargin{:});
component=p.Results.component;
timeincs=p.Results.timeincs;        % To be reported
timeinc=zeros(1,length(timesteps)); % Actually obtained
basecoord=p.Results.basecoordinate;
% for loop within postprocess(timesteps)
figure
hold on
for istep=1:length(timesteps)
    % pass variable to postprocess.plot method
    post=obj.Postprocess(timesteps(istep));
    timeinc(istep)=post.Inc;
    if strcmp(variable,'crackgeo')
            deformedPlus=post.EnrichItems{crackid}.IntPoints+post.EnrichItems{crackid}.Uplus;
            deformedMinus=post.EnrichItems{crackid}.IntPoints+post.EnrichItems{crackid}.Uminus;
            % one curve
            X=[deformedPlus(:,1);flipud(deformedMinus(:,1))];
            Y=[deformedPlus(:,2);flipud(deformedMinus(:,2))];
%             Y=(Y-basecoord)*1000;
            plot(X,Y,'LineWidth',2);
%             % two lines
%             plot(deformedPlus(:,1),deformedPlus(:,2),'LineWidth',1);
%             hold on
%             plot(deformedMinus(:,1),deformedMinus(:,2),'LineWidth',1);
            % fill plot
%             crackfill=[deformedPlus;flipud(deformedMinus)];
%             fill(crackfill(:,1),crackfill(:,2),'w','LineWidth',1)     % the mesh is gray, the crack space is white    
    else
        post.plot3d(variable,component);      
    end
    hold on
end
% add illustrations
titlestr=strcat('Shape Evolution of crack-', num2str(crackid));
title(titlestr)
xlabel('Distance from crack mouth (m)')
ylabel('Crack lip displacement (mm)')
%titlestr=['Snapshots of the ', variable,' at inquired time in the unit of second'];
%title(titlestr);
% timeinc=obj.NewtonRaphson.Timeinc(timesteps);
if ~isempty(timeincs)
    timeinc=timeincs;
end
legend(num2str(timeinc(:)));
hold off

