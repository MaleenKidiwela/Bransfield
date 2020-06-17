% File to change Water Picker information to Sting Ray formatting.
clear all
clc

%% Load OBS locations and depths
name = srStation.name;
for i=1:length(name)
    n= name(i);
    n = n{1};
    nameint(i) = str2num(n(4:5));
end
obs.name = nameint';
obs.lat = srStation.latitude;
obs.lon = srStation.longitude;
obs.depth = -1*srStation.elevation;
clear n nameint name

%% load the Shot lines
shots = readtable('orca_tomo_shotfile_final.txt')

shot(13,29) = struct( 'number', [], 'datetime', [], 'lat',[], 'lon',[], ...
    'shiplat',[], 'shiplon',[], 'depth',[], 'linename',[]);

for i = [13:16,18:26]
    filename = ['tlPick_BRA',num2str(i), '_Pw.mat'];
    loaded_file = load(filename);

    shot(i).number     = shots.x_shotnumber;
    for ii = 1: length(shots.date)
        strin(ii)   = strcat(string(shots.date(i))," ",string(shots.time(i)));
        strr(ii)   = datetime("2019-01-26 20:59:04.896785",'InputFormat','yyyy-mm-dd hh:mm:ss');
    end 
    
    shot(i).datetime   = strin;
    shot(i).lat        = str2num(c(:,36:45));
    shot(i).lon        = str2num(c(:,47:58));
    shot(i).shiplat    = str2num(c(:,60:68));
    shot(i).shiplon    = str2num(c(:,70:80));
    shot(i).depth      = str2num(c(:,82:88));
    shot(i).linename   = c(:,89:99);
    
end

%% Create tlArrival.
% Combine all the Arrival Time files.
file_name = NaN(1,length(obs.name));

all_tlArrival.eventid = [];
all_tlArrival.station = [];
all_tlArrival.time    = [];
all_tlArrival.error   = [];
tlArrival.eventid     = [];
tlArrival.station     = [];
tlArrival.time        = [];
tlArrival.error       = [];

for j = 1:length(obs.name);  % All OBSs.
    file_name = ['/Volumes/research/users/awells3/OBS_&_Event_Locations/Water_Wave_Picks/waterPick_OBS_', int2str(j), '.mat'];
    %file_name = ['/Network/Servers/thielsen.uoregon.edu/Volumes/home/awells3/Desktop/OBS Location/waterPick_OBS_', int2str(j), '.mat'];
    load(file_name);
    jtlArrival.station      = NaN(length(waterPick(j).shot),1);
    jtlArrival.station(:,1) = waterPick(j).obsname;
    all_tlArrival.time    = [tlArrival.time;    waterPick(j).artime];    % Creates all of the arrival times (secs)
    all_tlArrival.eventid = [tlArrival.eventid; waterPick(j).shot]; % Creates all of the shot numbersind = find(~isnan(all_tlArrival.time(:,1)));    
    all_tlArrival.station = [tlArrival.station; jtlArrival.station]; 
    all_tlArrival.error   = [tlArrival.error;   waterPick(j).sdv];   % Creates all of the arrival times (secs)
    ind = find(~isnan(all_tlArrival.time));
    tlArrival.time    = all_tlArrival.time(ind, 1);   % gets rid of NaNs in Matrices
    tlArrival.eventid = all_tlArrival.eventid(ind, 1); % gets rid of NaNs in Matrices
    tlArrival.station = all_tlArrival.station(ind, 1);
    tlArrival.error   = all_tlArrival.error(ind, 1); % gets rid of NaNs in Matrices
    
end

ind65= find(tlArrival.station==65);   % Find indices of OBS 15A.
ind66= find(tlArrival.station==66);   % Find indices of OBS 15B.
ind67= find(tlArrival.station==67);   % Find indices of OBS 15C.
ind68= find(tlArrival.station==68);   % Find indices of OBS 46A.

tlArrival.station    = int2str(tlArrival.station);
tlArrival.station    = cellstr(tlArrival.station);	% character (cell array)
tlArrival.station(ind65,1) = cellstr('15A');
tlArrival.station(ind66,1) = cellstr('15B');
tlArrival.station(ind67,1) = cellstr('15C');
tlArrival.station(ind68,1) = cellstr('46A');
tlArrival.station = strtrim(tlArrival.station); % Remove spaces
tlArrival.type       = ones(length(tlArrival.time),1);  % integer defining event type (Airguns)
tlArrival.phase      = cellstr(int2str(NaN(length(tlArrival.time),1)));  % Phase is 'Pw'.
tlArrival.phase(:,1) = cellstr('Pw');

cd '/Volumes/research/users/awells3/OBS_&_Event_Locations/OBS_&_Event_Relocation_Inversion/Inversion_Codes'
save('tlArrival.mat', 'tlArrival')  % Saving the results 


















