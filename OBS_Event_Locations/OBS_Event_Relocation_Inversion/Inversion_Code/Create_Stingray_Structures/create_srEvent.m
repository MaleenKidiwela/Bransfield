% File to change shot locations to Sting Ray formatting.
clear all
clc
addpath '/Volumes/research/users/awells3/OBS_&_Event_Locations/Water_Wave_Arrival_Picking_Codes'

%% Load OBS locations and depths
obs = read_obs_loc;

%% load the Shot lines
shot(1,64) = struct( 'number', [], 'datetime', [], 'lat',[], 'lon',[], ...
    'shiplat',[], 'shiplon',[], 'depth',[], 'linename',[]);
for i = [1:35,37:46];
    if i < 10;
        obsnum = ['0', int2str(i)];
    end
    if i >= 10;
        obsnum = int2str(i);
    end
    if i == 41;
        obsnum = '02A';
    end
    if i == 42;
        obsnum = '02R';
    end
    if i == 43;
        obsnum = '03R';
    end
    if i == 44;
        obsnum = '05A';
    end
    if i == 45;
        obsnum = '10A';
    end
    if i == 46;
        obsnum = '23R';
    end
    filename = ['/research/data/ETOMO/Data/obsip_shotlogs/MGL0910_', obsnum, '.obsip'];
    shot(i) = read_obsip_shotlog(filename);
end

%% Create srEvent.
srEvent.id        = [];
srEvent.latitude  = [];
srEvent.longitude = [];
srEvent.ot = [];

for iLine = 1:length(shot)
    srEvent.id         = [srEvent.id; shot(iLine).number]; % Create Shot Numbers
    srEvent.latitude   = [srEvent.latitude; shot(iLine).lat]; % Latitude of each Shot
    srEvent.longitude  = [srEvent.longitude; shot(iLine).lon];   % Longitude of each shot
    srEvent.ot = [srEvent.ot; shot(iLine).datetime]
end

srEvent.type = ones(length(srEvent.id),1);  % Type 1 signifies that airguns were used.

% Determine Airgun Depths.  

srEvent.elevation = NaN(length(srEvent.longitude),1);
depth1 = -0.015; % Depth airguns were dragged at for sequences 1-8. (km)
depth2 = -0.009; % Depth airguns were dragged for remaining sequences. (km)

for k = 1:length(srEvent.longitude);
    srEvent.origintime(k,1) = str2epoch(srEvent.ot(k,:));
    if srEvent.id(k) < 9000;
       srEvent.elevation(k,1) = depth1;
    else
       srEvent.elevation(k,1) = depth2;
    end
end
srEvent = rmfield(srEvent,'ot');

cd '/Volumes/research/users/awells3/OBS_&_Event_Locations/OBS_&_Event_Relocation_Inversion/Inversion_Codes'
save('srEvent.mat', 'srEvent')  % Saving the results 