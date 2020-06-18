function [x,y] = map2xy(mapx,mapy,srGeometry)

% MAP2XY -- Convert lat/lon positions to local cartesian
%           (Stingray utility)
%
%  [x,y] = map2xy(mapx,mapy,srGeometry)
%
%  Returns cartesian values in kilometers relative to the geographic center
%  of the experiment found in srGeometry.  Cartesian reference frame can be
%  rotated with respect to latitude and longitude or with respect to
%  easting and northing.  
%
%  INPUT:
%            mapx: either decimal longitude or easting
%            mapy: either decimal latitude or northing
%
%  OUTPUT:
%               x: x value in local cartesian, km
%               y: y value in local cartesian, km
%
%  srGeometry must have following fields for srGeometry.tf_latlon=true:
%
%       longitude: decimal degrees
%        latitude: decimal degrees
%        rotation: degrees
%       ellipsoid: radius and flattening from almanac
%
%  srGeometry must have following fields for srGeometry.tf_latlon=false:
%
%         easting: km
%        northing: km
%        rotation: degrees
%
%  Cartesian reference frame is right-handed, z positive upward, +'ve x is
%  to east if rotation is zer.  Length units of x and y are same as the
%  semimajor axis which is set in load_srGeometry using the alamanac (km).
%

%  Copyright 2010 Blue Tech Seismics, Inc.

%  Conversion and transformation done using matlab toolbox. Requires two
%  calls.  The geodetic292a toolbox by Michael Craymer has a function that
%  may do it in one call, if speed becomes an issue.

if srGeometry.tf_latlon == 1
    
    alt = 0;
    
    % Convert geodetic (ellipitcal) to ECEF
    
    [x,y,z] = geodetic2ecef(degtorad(mapy),degtorad(mapx),alt,srGeometry.ellipsoid);
    
    % Transform ECEF to ENU.  Origin of ENU coordinate system is
    % (srGeometry.latitude, srGeometry.longitude, 0)
    
    lat0rad   = degtorad(srGeometry.latitude);
    lon0rad   = degtorad(srGeometry.longitude);
    [dx,dy,~] = ecef2lv(x, y, z, lat0rad, lon0rad, alt, srGeometry.ellipsoid);
    
elseif srGeometry.tf_latlon== 0
    
    dx    = (mapx-srGeometry.easting);
    dy    = (mapy-srGeometry.northing );
    
end


rota  = srGeometry.rotation*pi/180;
sinrota = sin(rota);
cosrota = cos(rota);

x =  cosrota*dx + sinrota*dy;
y = -sinrota*dx + cosrota*dy;

