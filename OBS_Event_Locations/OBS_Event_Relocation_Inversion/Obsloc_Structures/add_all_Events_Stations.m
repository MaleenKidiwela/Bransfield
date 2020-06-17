%% File to add all events to Structure

%% Load all Stations & Events
addpath '/Volumes/research/users/awells3/OBS_&_Event_Locations/OBS_&_Event_Relocation_Inversion/Inversion_Codes'
srGeometry  = load_srGeometry('srGeometry');
srStation   = load_srStation('srStation',srGeometry);
srEvent     = load_srEvent('srEvent',srGeometry);
intsrStation = srStation;
intsrEvent = srEvent;

%% Load Structure to find station xyz and event xy errors.
load('Obsloc_structure_stations_xyz_events_xy');
a = s;
clear s;

%% Load Structure
load('Obsloc_structure_stations_xy_events_xy')

%% Change path to Obsloc_Structure
cd '/Volumes/research/users/awells3/OBS_&_Event_Locations/OBS_&_Event_Relocation_Inversion/Obsloc_Structures'

%% Find station xyz and event xy errors.
xStationerror = median(abs(a.xStation - s.xStation));
yStationerror = median(abs(a.yStation - s.yStation));
zStationerror = median(abs(a.zStation - s.zStation));
xEventerror   = median(abs(a.xEvent - s.xEvent));
yEventerror   = median(abs(a.yEvent - s.yEvent));
s.error_description = 'The errors for the stations and events were estimated from the median difference between the relocated station(xy) and event(xy) compared to relocated station(xyz) & events(xy).';

%% Plot initial shot locations
figure(1)
plot(srEvent.x,srEvent.y,'.');
hold on;

%% Change srEvents to include relocations.
ind = [];
srEvent.tf_relocate = zeros(srEvent.nevt,1);
for i = 1:srEvent.nevt
    if find(srEvent.id(i) == s.srEvent.id);
        ind = find(srEvent.id(i) == s.srEvent.id);
        srEvent.latitude(i) = s.srEvent.latitude(ind);
        srEvent.longitude(i) = s.srEvent.longitude(ind);
        srEvent.x(i) = s.srEvent.x(ind);
        srEvent.y(i) = s.srEvent.y(ind);
        srEvent.tf_relocate(i) = 1;
    end
end

% Plot all shots with new relocations
plot(srEvent.x,srEvent.y,'r.')
ylabel('Y, km');
xlabel('X, km');
title('All Event Locations: Initial vs Relocated (b. Inital, r. Relocated)');

% Change filename
srEvent.filename = 'srEvent_ETOMO_relocated';
srEvent.reloc_description = s.reloc_description;
srEvent.error_description = 'The errors for the events were estimated from the median of the difference between the relocated station(xy) and event(xy) compared to station(xyz) & events(xy)';
srEvent.xEventerror = xEventerror;
srEvent.yEventerror = yEventerror;
save srEvent_ETOMO_relocated srEvent

%% Plot initial shot locations
figure(2)
plot(srStation.x,srStation.y,'.');
hold on;

%% Change srStations to include relocations
% Missing instruments 13, 15, & 41.
srStation.tf_relocate = zeros(srStation.nsta,1);
for i = 1:srStation.nsta
    if i~=13 && i~=15 && i~=41;
        if i<13
            j = i;
        end
        if i ==14
            j = i-1;
        end
        if i > 14 && i<41
            j = i-2;
        end
        if i>41
            j = i-3;
        end
        srStation.latitude(i) = s.srStation.latitude(j);
        srStation.longitude(i) = s.srStation.longitude(j);
        srStation.elevation(i) = s.srStation.elevation(j);
        srStation.x(i) = s.srStation.x(j);
        srStation.y(i) = s.srStation.y(j);
        srStation.tf_relocate(i) = 1;
    end
end

% Plot all Stations with new relocations
plot(srStation.x,srStation.y,'r.')
ylabel('Y, km');
xlabel('X, km');
title('All Stations Locations: Initial vs Relocated (b. Inital, r. Relocated)');

% Save srStation
srStation.filename = 'srStation_ETOMO_relocated';
srStation.reloc_description = s.reloc_description;
srStation.error_description = 'The errors for the stations were estimated from the median of the difference between the relocated station(xy) and event(xy) compared to station(xyz) & events(xy)';
srStation.xStationerror = xStationerror;
srStation.yStationerror = yStationerror;
srStation.zStationerror = zStationerror;
save srStation_ETOMO_relocated srStation

%% Save final Obsloc structure
% Add all events and stations
s.xStationerror = xStationerror;
s.yStationerror = yStationerror;
s.zStationerror = zStationerror;
s.xEventerror   = xEventerror;
s.yEventerror   = yEventerror;
save Obsloc_structure s
    