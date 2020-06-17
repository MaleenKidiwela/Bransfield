function arrival = enter_arrival_time(file_name)

% Function to read the arrival time mat files

% INPUTS - the name of the arrival time file to read e.g 'waterPick_OBS_1_13-Jan-2010.mat'

% OUTPUTS
% Arrival time structure with the fields
%   obsname     - the OBS number
%   shot        - the shot number
%   artime      - the artime for a particular shot and OBS
%   sdv         - the standard deviation for artime

% Enter Arrival Time fiels


fid = fopen(file_name,'r');


header1 = fgetl(fid);        %#obsipshotfile v1.0 
header2 = fgetl(fid);        %#shotnumber date time sourceLat sourceLon shipLat shipLon waterDepth sciTag

c = fread(fid,'uchar');
c = setstr(c');
if length(file_name) == 68;
c = reshape(c,101,length(c)/101)';
end
if length(filename) == 57;
c = reshape(c,100,length(c)/100)';
end

shot.number     = str2num(c(:,1:7));
shot.datetime   = c(:,8:34);
shot.lat        = str2num(c(:,36:45));
shot.lon        = str2num(c(:,47:58));
shot.shiplat    = str2num(c(:,60:68));
shot.shiplon    = str2num(c(:,70:80));
shot.depth      = str2num(c(:,82:88));
shot.linename   = c(:,89:99);

fclose(fid);