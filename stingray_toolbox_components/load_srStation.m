function srStation=load_srStation(theFile,srGeometry)
% LOAD_srSTATTION -- Load the Stingray structure srStation
%                    (Stingray utility)
%
%  srStation = load_srStation(theFile,srGeometry)
%
%  Loads the srStation structure; checks required fields; sets derived
%  fields.
%
%  INPUT:
%               theFile:    Filename
%            srGeometry:    Stingray structure
%
%  OUTPUT:
%             srStation:    Stingray structure

%  Copyright 2010 Blue Tech Seismics, Inc.

load(theFile)

srStation.filename = theFile;

%%  Required fields

if isempty(srStation.name) || isempty(srStation.elevation) 
    error('Required fields missing in srStation')
end

%%  Derived fields

%  Check geographic vs. cartesian
if srGeometry.tf_latlon
    if isempty(srStation.latitude) || isempty(srStation.longitude)
        error('Required fields missing in srStation')
    end
    [srStation.x srStation.y] = map2xy(srStation.longitude, srStation.latitude, ...
        srGeometry);
elseif ~srGeometry.tf_latlon
    if isempty(srStation.easting) || isempty(srStation.northing)
        error('Required fields missing in srStation')
    end
    
    [srStation.x srStation.y] = map2xy(srStation.easting, srStation.northing, ...
        srGeometry);
elseif isempty(srGeometry.tf_latlon)
    error('ERROR:  srGeometry.tf_latlon.tf_laton has no value')
end

srStation.nsta = length(srStation.name);

display(srStation);


