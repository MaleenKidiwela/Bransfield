%% Compare Different Structures.

% Load Structures to Compare
addpath('/Volumes/research/users/awells3/OBS_&_Event_Locations/OBS_&_Event_Relocation_Inversion/Obsloc_Structures');
load('Obsloc_structure_stations_xy_events_fixed')
a = s;
load('Obsloc_structure_stations_xyz_events_fixed')
b = s;
load('Obsloc_structure_stations_xyz_events_xy')
c = s;
load('Obsloc_structure_stations_xyz_events_xy_smallzerror')
d = s;
load('Obsloc_structure_stations_xy_events_xy');
e = s;
load('Obsloc_structure_stations_xy_events_xy_2increase');
f = s;
load('Obsloc_structure_stations_xy_events_xy_changedOBS45');
g = s;
load('Obsloc_structure_stations_xy_events_xy_2increase_changedOBS45');
h = s;

% Number of Stations and shots and indexing vectors
m = length(c.ievt);
n = length(c.ista);

aMatrix = a.time - a.timePred;
bMatrix = b.time - b.timePred;
cMatrix = c.time - c.timePred;
dMatrix = d.time - d.timePred;
eMatrix = e.time - e.timePred;
fMatrix = f.time - f.timePred;
gMatrix = g.time - g.timePred;
hMatrix = h.time - h.timePred;

%% Plot changes for the final iteration for changes in z and no changes in z
% % 
% figure(1); clf
% plot(a.ista,(abs(a.xStation_change)+abs(a.yStation_change)+abs(a.zStation_change))/3,'*b')
% hold on
% % plot(b.ista,(abs(b.xStation_change)+abs(b.yStation_change))/2,'*g')
% grid on
% title('Average Change in Station location from initial locations (Change in z - blue, No Change in z - green)')
% xlabel('Station Index')
% ylabel('km');

% figure(2);
% matrix = [((abs(a.xStation_change)+abs(a.yStation_change)+abs(a.zStation_change))/3), ((abs(b.xStation_change)+abs(b.yStation_change))/2)];
% subplot(211)
% hist(matrix);
% title('Average Change in Station location from initial locations (Change in z - blue, No Change in z - red)')
% xlabel('Ave Change in Station Location - km')
% ylabel('Number of Stations');
% subplot(212)
% plot(matrix,'*');
% title('Average Change in Station location from initial locations (Change in z - blue, No Change in z - green)')
% xlabel('Station Number')
% ylabel('Change in Location - km');

% figure(3);
% subplot(211)
% hist(abs(aMatrix));
% title('Travel time misfits (Observed minus Predicted) (Change in z)')
% ylabel('Event Index')
% xlabel('Misfit, s');
% subplot(212)
% % hist(abs(bMatrix));
% title('Travel time misfits (Observed minus Predicted) (No Change in z)')
% ylabel('Event Index')
% xlabel('Misfit, s');
% hist(abs(aMatrix));


for j = 1:length(a.ista);
    stdaMatrix(j) = std(abs(aMatrix(~isnan(aMatrix(:,j)),j)));
    stdbMatrix(j) = std(abs(bMatrix(~isnan(bMatrix(:,j)),j)));
    stdcMatrix(j) = std(abs(cMatrix(~isnan(cMatrix(:,j)),j)));
    stddMatrix(j) = std(abs(dMatrix(~isnan(dMatrix(:,j)),j)));
    stdeMatrix(j) = std(abs(eMatrix(~isnan(eMatrix(:,j)),j)));
    stdfMatrix(j) = std(abs(fMatrix(~isnan(fMatrix(:,j)),j)));
    stdgMatrix(j) = std(abs(gMatrix(~isnan(gMatrix(:,j)),j)));
    meanaMatrix(j) = mean(abs(aMatrix(~isnan(aMatrix(:,j)),j)));
    meanbMatrix(j) = mean(abs(bMatrix(~isnan(bMatrix(:,j)),j)));
    meancMatrix(j) = mean(abs(cMatrix(~isnan(cMatrix(:,j)),j)));
    meandMatrix(j) = mean(abs(dMatrix(~isnan(dMatrix(:,j)),j)));
    meaneMatrix(j) = mean(abs(eMatrix(~isnan(eMatrix(:,j)),j)));    
    meanfMatrix(j) = mean(abs(fMatrix(~isnan(fMatrix(:,j)),j)));    
    meangMatrix(j) = mean(abs(gMatrix(~isnan(gMatrix(:,j)),j)));
    meanhMatrix(j) = mean(abs(hMatrix(~isnan(hMatrix(:,j)),j)));
end

%% Compare Means
[k, ind] = sort(meaneMatrix)
figure(10)
plot(meaneMatrix(ind),'b.-')
hold on
plot(meanfMatrix(ind),'r.-')
plot(meangMatrix(ind),'g.-')
plot(meanhMatrix(ind),'m.-')
plot(meancMatrix(ind),'c.-')
legend('stations xy events xy','stations xy events xy 2increase','stations xy events xy changedOBS45','stations xy events xy 2increase changedOBS45','stations xyz events xy')
xlabel('stations')
ylabel('mean misfit (time-timePred) (sec)')
title('Mean Misfit Comparison between Structures')
hold off

%% Plot mean travel time misfits per station on  map
figure(1)
edata = meaneMatrix;%(:,ind);
[o,n]=size(edata); 
edata = reshape(edata,o*n,1);
edata = edata(find(~isnan(edata)));
lim = [min(edata) max(edata)];
col = LinColor(edata);
for i=1:length(e.ista)%[1:41,43:50,52:65];%1:length(data);
plot(e.xStation(i),e.yStation(i),'o','MarkerFaceColor',col(i,:),'MarkerEdgeColor',col(i,:),'Markersize',8);
hold on;
text(e.xStation(i)+1,e.yStation(i)+1,cell2mat(e.srStation.name(i)),'fontsize',10)
text(e.xStation(i)+1,e.yStation(i)-.5,num2str(meaneMatrix(i),2))
end
axis('equal');
ylabel('Y, km');
xlabel('X, km');
title(['Mean Misfit for all Stations (Obsloc structure stations xy events xy) (seconds)']);
caxis(lim);
hcb = colorbar;
h = gca;
set(gcf,'currentaxes',hcb);
title('Residual, s');
set(gcf,'currentaxes',h);

figure(2)
bdata = meanbMatrix;
[o,n]=size(bdata); 
bdata = reshape(bdata,o*n,1);
bdata = bdata(find(~isnan(bdata)));
lim = [min(bdata) max(bdata)];
col = LinColor(bdata);

for i=1:length(b.ista)
plot(b.xStation(i),b.yStation(i),'o','MarkerFaceColor',col(i,:),'MarkerEdgeColor',col(i,:),'Markersize',8);
hold on;
text(b.xStation(i)+1,b.yStation(i)+1,cell2mat(b.srStation.name(i)),'fontsize',10)
text(b.xStation(i)+1,b.yStation(i)-.5,num2str(meanbMatrix(i),2))
end
axis('equal');
ylabel('Y, km');
xlabel('X, km');
title(['Mean Misfit for all Stations (Obsloc structure stations xyz events fixed 5increase) (seconds)']);
caxis(lim);
hcb = colorbar;
h = gca;
set(gcf,'currentaxes',hcb);
title('Residual, s');
set(gcf,'currentaxes',h);


figure(3)
cdata = meancMatrix;
[o,n]=size(cdata); 
cdata = reshape(cdata,o*n,1);
cdata = cdata(find(~isnan(cdata)));
lim = [min(cdata) max(cdata)];
col = LinColor(cdata);
for i=1:length(c.ista)
plot(c.xStation(i),c.yStation(i),'o','MarkerFaceColor',col(i,:),'MarkerEdgeColor',col(i,:),'Markersize',8);
hold on;
text(c.xStation(i)+1,c.yStation(i)+1,cell2mat(c.srStation.name(i)),'fontsize',10)
text(c.xStation(i)+1,c.yStation(i)-.5,num2str(meancMatrix(i),2))
end
axis('equal');
ylabel('Y, km');
xlabel('X, km');
title(['Mean Misfit for all Stations (Obsloc structure stations xyz events xy) (seconds)']);
caxis(lim);
hcb = colorbar;
h = gca;
set(gcf,'currentaxes',hcb);
title('Residual, s');
set(gcf,'currentaxes',h);

% figure(7)
% data = stdaMatrix;%(:,ind);
% [m,n]=size(data); 
% data = reshape(data,m*n,1);
% data = data(find(~isnan(data)));
% lim = [min(data) max(data)];
% col = LinColor(data);
% for i=1:length(a.ista)%[1:41,43:50,52:65];%1:length(data);
% plot(a.xStation(i),a.yStation(i),'o','MarkerFaceColor',col(i,:),'MarkerEdgeColor',col(i,:),'Markersize',8);
% hold on;
% text(a.xStation(i)+1,a.yStation(i)+1,cell2mat(a.srStation.name(i)))
% text(a.xStation(i)+1,a.yStation(i)-.5,num2str(stdaMatrix(i),2))%, ...
% end
% axis('equal');
% ylabel('Y, km');
% xlabel('X, km');
% title(['STD for Station Misfits']);
% caxis(lim);
% hcb = colorbar;
% h = gca;
% set(gcf,'currentaxes',hcb);
% title('STD');
% set(gcf,'currentaxes',h);

%% Plot comparison of Structures
% Mean Misfits
figure(5)
[B,ind1] =sort(meancMatrix);
plot(meanaMatrix(ind1),'k.')
hold on
plot(meanbMatrix(ind1),'r.')
plot(meancMatrix(ind1),'g.')
plot(meandMatrix(ind1),'b.')
plot(meaneMatrix(ind1),'m.')
% plot(meanfMatrix(ind1),'c.')
% plot(meangMatrix(ind1),'y.')
ylabel('mean misfit (s)');
xlabel('stations');
title(['Mean Misfit for all Stations (seconds)']);
hold off
legend('stations xy events fixed','stations xyz events fixed','stations xyz events xy','stations xyz events xy zErrorApr = 0.001','stations xy events xy')
%,'stations xyz events fixed 5increase','stations xyz events xy 8increase','stations xyz events fixed 8increase')
grid on

% Station Depths
[B,ind] =sort(e.zStation);
figure(6)
plot(e.zStation(ind),'b*','MarkerSize',10)
hold on
plot(f.zStation(ind),'r.')
plot(g.zStation(ind),'g.')
plot(h.zStation(ind),'m.')
plot(c.zStation(ind),'c.')
% plot(f.zStation(ind),'c.')
% plot(g.zStation(ind),'y.')
ylabel('depths (km)');
xlabel('stations');
title(['Station Depths (seconds)']);
hold off
legend('stations xy events xy','stations xy events xy 2increase','stations xy events xy changedOBS45','stations xy events xy 2increase changedOBS45','stations xyz events xy')
grid on
set(gca,'YDir','reverse')

% Station Depth Differences
figure(7)
%plot(a.zStation(ind),'k*','MarkerSize',10)
plot(abs(b.zStation-a.zStation),'r.-')
hold on
plot(abs(c.zStation-a.zStation),'g.-')
% plot(abs(d.zStation-a.zStation),'b.-')
% plot(abs(e.zStation-a.zStation),'m.-')
% plot(abs(f.zStation-a.zStation),'c.-')
% plot(abs(g.zStation-a.zStation),'y.-')
ylabel('depths difference (km)');
xlabel('stations');
title(['Station Depths Differences (seconds)']);
hold off
legend('stations xyz events fixed','stations xyz events xy','stations xyz events xy 5increase','stations xyz events fixed 5increase','stations xyz events xy 8increase','stations xyz events fixed 8increase')
grid on

% Plot Map of Relocated Shots & Stations compared to Original
figure(8)
plot(e.xStation_int,e.yStation_int,'rd','MarkerSize',10)
hold on
plot(c.xStation,c.yStation,'gd','MarkerSize',10)
plot(e.xStation,e.yStation,'bd','MarkerSize',10)
plot(e.xEvent_int,e.yEvent_int,'r.-')
plot(c.xEvent,c.yEvent,'g.-')
plot(e.xEvent,e.yEvent,'b.-')
ylabel('Y (km)');
xlabel('X (km)');
title(['Map of Station and Event Relocations (r inital loc, b stations xy events xy, g stations xyz events xy)']);
hold off


%% Check Time vs Distance
% [shotIndex,StationIndex] = find(ones(m,n));
% DXa = a.srEvent.x(shotIndex,:) - a.srStation.x(StationIndex,:);
% DYa = a.srEvent.y(shotIndex,:) - a.srStation.y(StationIndex,:);
% rangea = sqrt( DXa.^2 + DYa.^2 );
% deptha = a.srStation.z(StationIndex,:);
% gundeptha = a.srEvent.elevation(shotIndex,:);
% da = sqrt(rangea.^2 + (deptha+gundeptha).^2) ;  % the distance
% da = reshape( da ,m,n );  % the distance
% rangea = reshape( rangea ,m,n );  % the distance
% figure(8)
% plot(rangea,a.time,'.'), hold on
% DXc = c.srEvent.x(shotIndex,:) - c.srStation.x(StationIndex,:);
% DYc = c.srEvent.y(shotIndex,:) - c.srStation.y(StationIndex,:);
% rangec = sqrt( DXc.^2 + DYc.^2 );
% depthc = c.srStation.z(StationIndex,:);
% gundepthc = c.srEvent.elevation(shotIndex,:);
% dc = sqrt(rangec.^2 + (depthc+gundepthc).^2) ;  % the distance
% dc = reshape( dc ,m,n );  % the distance
% rangec = reshape( rangec,m,n );  % the distance
% figure(8)
% plot(rangec,c.time-rangec/1.5,'.'), hold on
% xlim([0 12])
% xlabel('Range (km)');
% ylabel('Time-(Range/1.5) (s)');
% title(['Time vs Range']);
% text(6,1.2, '. stations xyz events xy');
% text(6,1.1, 'o stations xyz events xy 5increase');
% text (6,1, '* stations xyz events xy 8increase');
% grid on

% DXd = d.srEvent.x(shotIndex,:) - d.srStation.x(StationIndex,:);
% DYd = d.srEvent.y(shotIndex,:) - d.srStation.y(StationIndex,:);
% ranged = sqrt( DXd.^2 + DYd.^2 );
% depthd = d.srStation.z(StationIndex,:);
% gundepthd = d.srEvent.elevation(shotIndex,:);
% dd = sqrt(ranged.^2 + (depthd+gundepthd).^2) ;  % the distance
% dd = reshape( dd ,m,n );  % the distance
% ranged = reshape( ranged ,m,n );  % the distance
% plot(ranged,d.time-ranged/1.5,'o');
% 
% DXf = f.srEvent.x(shotIndex,:) - f.srStation.x(StationIndex,:);
% DYf = f.srEvent.y(shotIndex,:) - f.srStation.y(StationIndex,:);
% rangef = sqrt( DXf.^2 + DYf.^2 );
% depthf = d.srStation.z(StationIndex,:);
% gundepthf = f.srEvent.elevation(shotIndex,:);
% df = sqrt(rangef.^2 + (depthf+gundepthf).^2) ;  % the distance
% df = reshape( df ,m,n );  % the distance
% rangef = reshape( rangef ,m,n );  % the distance
% plot(rangef,f.time-rangef/1.5,'*');
% hold off

% figure(9)
% plot(c.xStation_int,c.yStation_int,'bd','MarkerSize',10), hold on
% plot(c.xStation,c.yStation,'rs','MarkerSize',10)
% plot(d.xStation,d.yStation,'gs','MarkerSize',10)
% plot(f.xStation,f.yStation,'ms','MarkerSize',10)
% %plot(a.xStation,a.yStation,'cs','MarkerSize',10)
% plot(c.xEvent_int,c.yEvent_int,'b*-')
% plot(c.xEvent,c.yEvent,'r.-')
% plot(d.xEvent,d.yEvent,'g.-')
% plot(f.xEvent,f.yEvent,'m.-')
% xlabel('X (km)');
% ylabel('Y (km)');
% title(['Station and Event Locations (Reloc stations xyz events xy) (b initial loc, r original reloc, g 5increase , m 8increase)']);
% hold off
