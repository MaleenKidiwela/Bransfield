% File to change OBS locations to Sting Ray formatting.
clear all
clc

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


%% Create srStation.  
% OBS S15A is obs.name 65
% OBS S15B is obs.name 66
% OBS S15C is obs.name 67
% OBS S46A is obs.name 68
for k = 1:64;
    srStation.name(k) = cellstr(int2str(obs.name(k)));
end

srStation.name(65)  = cellstr('15A');
srStation.name(66)  = cellstr('15B');
srStation.name(67)  = cellstr('15C');
srStation.name(68)  = cellstr('46A');
srStation.name      = srStation.name';
srStation.longitude = obs.lon;
srStation.latitude  = obs.lat;
srStation.elevation = -0.001*obs.depth;  % Kilometers relative to sea level.  Positive values are above sea level.
srStation.filename  = ['/Network/Servers/thielsen.uoregon.edu/Volumes/home/awells3/Desktop/OBS Location/OBS_location.txt'];

cd '/Volumes/research/users/awells3/OBS_&_Event_Locations/OBS_&_Event_Relocation_Inversion/Inversion_Codes'
save('srStation.mat', 'srStation')  % Saving the results 