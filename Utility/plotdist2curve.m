plot(tempsegment(:,1),tempsegment(:,2),'k-o',plist(:,1),plist(:,2),'r*')
hold on
plot(xy(:,1),xy(:,2),'g*')
line([plist(:,1),xy(:,1)]',[plist(:,2),xy(:,2)]','color',[0 0 1])
axis equal