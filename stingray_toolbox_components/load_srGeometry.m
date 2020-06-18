function srGeometry=load_srGeometry(theFile)
% LOAD_srGEOMETRY -- Load the Stingray structure srGeometry
%                   (Stingray utility)
%
%  srGeometry = load_srGeometry(theFile)
%
%  Loads the srGeometry structure; checks required fields; sets derived
%  fields.
%
%  INPUT:
%               theFile:    Filename
%
%  OUTPUT:
%            srGeometry:    Stingray structure

%  Copyright 2010 Blue Tech Seismics, Inc.

load(theFile)
srGeometry.filename = theFile;

%% Required fields

if isempty(srGeometry.tf_latlon)
    error('Must set srGeometry.tf_latlon')
elseif srGeometry.tf_latlon
    if isempty(srGeometry.latitude) || isempty(srGeometry.longitude) || ...
            isempty(srGeometry.rotation)
        error('Required variable missing in srGeometry')
    end
elseif ~srGeometry.tf_latlon
     if isempty(srGeometry.easting) || isempty(srGeometry.northing) || ...
            isempty(srGeometry.rotation)
        error('Required variable missing in srGeometry')
     end 
end       

%% Derived fields

% Use the almanac to return ellipsoid values in kilometers

if srGeometry.tf_latlon
    srGeometry.ellipsoid = almanac('earth','wgs84','kilometers');
end

display(srGeometry)


