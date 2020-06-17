%% Check Sounding Velocities for the OBS Relocations

%% Load Sounding Velocity files
load /Volumes/research/users/awells3/XBT/SoundingVel_09m.dat
load /Volumes/research/users/awells3/XBT/SoundingVel_15m.dat
load /Volumes/research/users/awells3/XBT/AverageSoundingVel.dat   

% Load Obsloc Structure - Relocated Stations (xy), fixed Events
addpath('/Volumes/research/users/awells3/OBS_&_Event_Locations/OBS_&_Event_Relocation _Inversion/Obsloc_Structures')
load('Obsloc_structure_stations_xy_events_fixed')
%load('Obsloc_structure_stations_xy_events_xy');

% Set Range for shots to examine around each OBS.
max_range = 2; % Outside range - km
min_range = 1; % Inside range - km

%% Plot Average Sounding Velocities
% figure(1);
% plot(AverageSoundingVel(:,2),-AverageSoundingVel(:,1),'-m','linewidth',2);
% hold on
% plot(SoundingVel_09m(:,2),-SoundingVel_09m(:,1),'-c','linewidth',2); %Sounding Vel from 9m
% plot(SoundingVel_15m(:,2),-SoundingVel_15m(:,1),'-k','linewidth',2); %Sounding vel from 15m
% hold off

%% Load Velocity Tables from William Wilcock.
addpath('/Volumes/research/users/awells3/OBS_&_Event_Locations/OBS_&_Event_Relocation _Inversion/ray_tables')
p.nlookup = 2;
%p.lookup(1).file = 'ray_ttTable_AverageWaterVel1_9.mat'; 
p.lookup(1).file = 'ray_ttTable_AverageWaterVel1_5increase_9.mat'; 
%p.lookup(1).file = 'TTtables_simple_9m.mat';
p.lookup(1).eventElevation = -0.009;
%p.lookup(2).file = 'ray_ttTable_AverageWaterVel1_15.mat'; 
p.lookup(2).file = 'ray_ttTable_AverageWaterVel1_5increase_15.mat'; 
%p.lookup(2).file = 'TTtables_simple_15m.mat';
p.lookup(2).eventElevation = -0.015;

lookup = ones(length(s.ievt),1);
for i = 1:p.nlookup
  eval(['load ' p.lookup(i).file])
  p.lookup(i).ttTable   = ttTable;
  p.lookup(i).dtdrTable = dtdrTable;
  p.lookup(i).dtdzTable = dtdzTable;
  p.lookup(i).rTable    = rTable;
  p.lookup(i).zTable    = zTable;
  if i>1
    lookup(s.srEvent.elevation==p.lookup(i).eventElevation) = i;
  end
end

% Number of Stations and shots and indexing vectors
m = length(s.ievt);
n = length(s.ista);
[shotIndex,StationIndex] = find(ones(m,n));
DX = s.srEvent.x(shotIndex,:) - s.srStation.x(StationIndex,:);
DY = s.srEvent.y(shotIndex,:) - s.srStation.y(StationIndex,:);
range = sqrt( DX.^2 + DY.^2 );
depth = s.srStation.z(StationIndex,:);
d = sqrt(range.^2 + depth.^2) ;  % the distance
d = reshape( d ,m,n );  % the distance
timePred = NaN(m,n);
dtdr = NaN(m,n);
dtdz = NaN(m,n);
srEventid = s.srEvent.id(shotIndex,:);
srStationid = s.srStation.name(StationIndex,:);

for ilookup = 1:p.nlookup
    timePred1 = interp2(p.lookup(ilookup).zTable,p.lookup(ilookup).rTable,...
        p.lookup(ilookup).ttTable,depth,range,'*linear',NaN);
    timePred1 = reshape(timePred1,m,n);
    timePred(lookup==ilookup,:) = timePred1(lookup==ilookup,:);
    
    dtdr1 = interp2(p.lookup(ilookup).zTable,p.lookup(ilookup).rTable,...
        p.lookup(ilookup).dtdrTable,depth,range,'*linear',NaN);
    dtdr1 = reshape(dtdr1,m,n);
    dtdr(lookup==ilookup,:) = dtdr1(lookup==ilookup,:);
    
    dtdz1 = interp2(p.lookup(ilookup).zTable,p.lookup(ilookup).rTable,...
        p.lookup(ilookup).dtdzTable,depth,range,'*linear',NaN);
    dtdz1 = reshape(dtdz1,m,n);
    dtdz(lookup==ilookup,:) = dtdz1(lookup==ilookup,:);
end

% Calculate Depth from Airgun
for i =1:length(s.ista);
    d(:,i) = d(:,i) + s.srEvent.elevation;
end

% Calculate Sounding Velocities from Travel Time Tables
table_vel = d./timePred;

%% Calculate ETOMO Survey's Sounding Velocities for each OBS 

% Find Events within Specified Range to OBS.
Jin = find(range < max_range);
Jrange = range(Jin(find(range(Jin) > min_range)));
for i = 1:length(Jrange)
    J(i) = find(Jrange(i)==range); 
end

% Calculate Sounding Velocities for Specified Range
vel = nan(1,length(range));
vel(J) = (d(J)./s.time(J)); %km/s
vel(1,(find(vel(J(length(J)))==vel)+1):n*m) = NaN;
jj = J(~isnan(vel(J)))
median_vel = median(vel(jj));
median_depth = depth(find(vel == median_vel));
vel = reshape(vel,m,n);

% Calculate Average Sounding Velocities from Experiment and Travel Time Tables.
% Plot Comparison Graphs
for i = 1:n % For each OBS
    OBS_vel = vel(~isnan(vel(:,i)),i);
    ave_vel(i) = mean(vel(~isnan(vel(:,i)),i))';
    OBS_table_vel = table_vel(~isnan(vel(:,i)),i);
    ave_table_vel(i) = mean(table_vel(~isnan(vel(:,i)),i))';
    figure(2)
    plot(i*ones(size(OBS_vel)),OBS_vel.*1000,'g.')
    hold on
    plot(i*ones(size(OBS_vel)),OBS_table_vel.*1000,'y.')
    plot(i,ave_vel(i)*1000,'b.')
    plot(i,ave_table_vel(i)*1000,'r.','markersize',5)
    text(i-.3,ave_vel(i)+.01,cell2mat(s.srStation.name(i)),'fontsize',10)
end
title(['\fontsize{12}Sounding Velocity (', int2str(min_range),' to ', int2str(max_range),'km range) (OBS vel -g, Ave OBS vel -b, Table vel -y, Ave Table vel -r)'])
xlabel(['\fontsize{12}OBS'])
ylabel(['\fontsize{12}Velocity m/s']);
hold off


%% Plot Velocities to Compare
figure(7)
plot(SoundingVel_09m(:,2),SoundingVel_09m(:,1)./1000,'-m','linewidth',1.5); %Sounding Vel from 9m
hold on;
plot(SoundingVel_15m(:,2),SoundingVel_15m(:,1)./1000,'-g','linewidth',2); %Sounding vel from 15m
plot(ave_table_vel.*1000,s.srStation.z,'r.',ave_vel.*1000,s.srStation.z,'b.')
plot(median_vel*1000, median_depth, '*b', 'MarkerSize',10)
for i=1:length(s.ista);
text(ave_vel(i).*1000+.3,s.srStation.z(i)+.01,cell2mat(s.srStation.name(i)),'fontsize',10)
end
set(gca,'YDir','reverse')
title(['\fontsize{12}Sounding Velocity (', int2str(min_range),' to ', int2str(max_range), 'km range) (XBT 9m -m, XBT 15m -g, Calc SV b., Table Lookup r.)'])
xlabel(['\fontsize{12}Velocity m/s'])
ylabel(['\fontsize{12}Depth of OBSs - km']);
hold off

%% Plot mean velocities per station on  map
figure(4)
data = ave_vel;
[m,n]=size(data);
data = reshape(data,m*n,1);
id_nodata = find(isnan(data));
data = data(find(~isnan(data)));
lim = [min(data) max(data)];
col = LinColor(data);
for i=1:(n-length(id_nodata))
    j = i;
    if id_nodata;
        if i>=id_nodata;
            j = i + 1;
        end
    end
    plot(s.xStation(j),s.yStation(j),'o','MarkerFaceColor',col(i,:),'MarkerEdgeColor',col(i,:),'Markersize',8);
    hold on;
    text(s.xStation(j)+1,s.yStation(j)+1,cell2mat(s.srStation.name(j)),'fontsize',10)
    text(s.xStation(j)+1,s.yStation(j)-.5,num2str(ave_vel(j),3))
end

axis('equal');
ylabel('Y, km');
xlabel('X, km');
title(['\fontsize{12}Mean Sounding Velocities for all Stations (', int2str(min_range),' to ', int2str(max_range),'km range) (km/s)']);
caxis(lim);
hcb = colorbar;
h = gca;
set(gcf,'currentaxes',hcb);
title('\fontsize{12}Vel, km/s');
set(gcf,'currentaxes',h);
hold off

