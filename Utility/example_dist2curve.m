x=1:0.2:10;
y=sin(x);
curvexy=[x;y]';
xm=1:0.1:10;
ym=0.5*ones(1,length(xm));
mapxy=[xm;ym]';
plot(x,y,'-')
hold on
plot(xm,ym,'o')
hold off
[xy,dist,t,minind]=distance2curve(curvexy,mapxy,'linear');
