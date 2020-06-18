% Creates experiment basemap and plots loaded shots if requested 
figure(100)
clf
plot(srEvent.x,srEvent.y,'.k','markersize',4);
hold on
title('Experiment Geometry')
xlabel('X, km')
ylabel('Y, km')
legend('Event','Station','location','bestoutside')
axis('equal')
        
if currentMenu.geometryPlot.loaded    
  i = unique(traceMetaData.eventIndex);
  hGP1 = plot(srEvent.x(i),srEvent.y(i),'or','markersize',4,'markerfacecolor','r');
else
  clear hGP1 hGP2
end

plot(srStation.x,srStation.y,'sk','markerfacecolor','k','markersize',8);

if currentMenu.geometryPlot.loaded    
  i = unique(traceMetaData.stationIndex);
  hGP2 = plot(srStation.x(i),srStation.y(i),'sk','markerfacecolor','c','markersize',12);
  legend('Event','Loaded Event','Station','Loaded Station','location','bestoutside')
else
  clear hGP1 hGP2
end

clear hGP3 hGP4
figure(hFig)