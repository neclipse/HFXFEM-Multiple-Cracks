clear;clc;
x = [0 1 2 2 1 1 0 0];
y = [0 0 0 1 1 2 2 1];
for i = 1 : length(x)
    text(x(i)+0.1,y(i)+0.1,sprintf('%d',i))
    hold on
end
quad = [1 2 5 8;
    2 3 4 5;
    8 5 6 7];

plot(x,y,'k.-')
axis off

z = sin(pi*x).*cos(pi*y);
c = z;
figure

ax=newplot;
fc = get(gcf,'Color');
h = patch('faces',quad,'vertices',[x(:) y(:) z(:)],'facevertexcdata',c(:),...
    'facecolor',fc,'edgecolor',get(ax,'defaultsurfacefacecolor'),...
    'facelighting', 'none', 'edgelighting', 'flat',...
    'parent',ax);