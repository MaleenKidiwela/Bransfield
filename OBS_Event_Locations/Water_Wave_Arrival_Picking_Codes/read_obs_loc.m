 function obs = read_obs_loc(filename)

% Function to read the obs_loc file and return the data below header 
% in the structure shot with the fields.  
%
% INPUTS - the name of the obsip file to read e.g 'MGL0910_09.obsip'
%
% OUTPUTS
% OBS structure with the fields
%   NAME        - the OBS number
%   LAT         - the latitude of the shot
%   LON         - the longitude of the shot
%   DEPTH       - the water depth of the OBS


% Find OBS locations.
% OBS S15A is obs.name 65
% OBS S15B is obs.name 66
% OBS S15C is obs.name 67
% OBS S46A is obs.name 68


fid = fopen('OBS_location.txt  ','r');


nsta     = fscanf(fid,'%i',1); 
statData = fscanf(fid,'%i %g %g %g',[4,nsta]);
statData = statData';
    
obs.name        = statData(:,1);
obs.lat         = statData(:,2);
obs.lon         = statData(:,3);
obs.depth       = statData(:,4);

fclose(fid);

