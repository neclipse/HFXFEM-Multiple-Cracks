function showme( obj,typex,varargin )
% Show the most requested plots from CrackBody results with annotations
%   Input explanation
% typex: % 1-x coordinates of IntPoints; 2-y coordinates of IntPoints; 3- Real distance from one end
% first: first (primary variable)
p=inputParser;
% addRequired(p,'typex');
addRequired(p,'first');
addOptional(p,'second','Pfrack');
addParameter(p,'component',1);
parse(p,varargin{:});
% typex=p.Results.typex;
first=p.Results.first;
second=p.Results.second;
component=p.Results.component;
% determine xlist
switch typex
    case 1
        x=obj.IntPoints(:,1);
        xlab='X-coordinates of gaussian points along the crack (m)';
    case 2
        x=obj.IntPoints(:,2);
        xlab='Y-coordinates of gaussian points along the crack (m)';
    case 3
        tip=obj.Mygeo.Tips(2,:);
        % for now
        x=sqrt((obj.IntPoints(:,1)-tip(1)).^2+(obj.IntPoints(:,2)-tip(2)).^2);
        xlab='Linear distance between the gaussian points from left crack tip (m)';
end

% withdraw the first variable by eval
expression=strcat('obj.',first);
y1=eval(expression);
% Make plots
if ismember('second',p.UsingDefaults)
    % only plot the first variable
    figure
    plot(x,y1*1000,'LineWidth',2);
    xlabel(xlab,'FontSize',14);
    ylabel(first,'FontSize',14);
else
    % use plotyy to plot two variables
    exp2=strcat('obj.',second);
    y2temp=eval(exp2);
    y2=y2temp(:,component);
    figure
    [ax,h1,h2]=plotyy(x,y1*1000,x,y2*1000);
    set(ax,'FontSize',14)
    set(h1,'LineWidth',2);
    set(h2,'LineWidth',2);
    titlestr=strcat(second,' vs. ',first);
    title(titlestr,'FontSize',14);
    xlabel(xlab,'FontSize',14)
    yticks(ax(1),'auto')
    yticks(ax(2),'auto')
    ylabel(ax(1),first,'FontSize',14)
    ylabel(ax(2),second,'FontSize',14)
    ax(1).Position=[0.13,0.15,0.7250,0.7750];
%     ticks
end
end

