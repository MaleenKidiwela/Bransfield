function [mapx,mapy] = xy2map(xlc,ylc,srGeometry)

% XY2MAP -- Convert local cartesian to geodetic or UTM coordinates
%           (Stingray utility)
%
% [mapx,mapy] = xy2map(xlc,ylc,srGeometry)
%
%  Local cartesian reference frame can be rotated with respect to latitude
%  and longitude or with respect to easting and northing.  
%
%  INPUT:
%             xlc: x value in local cartesian, km
%             ylc: y value in local cartesian, km 
%      srGeometry: Stingray structure
%
%  OUTPUT:
%            mapx: either decimal longitude or easting
%            mapy: either decimal latitude or northing
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
%  to right.  Length units of x and y are same as the semimajor axis which
%  is set in load_srGeometry using the alamanac (km).
%

%  Copyright 2010 Blue Tech Seismics, Inc.

%  Transform from local cartesian to ENU.  Up value will be set to zero.

rad = pi/180.;
csr = cos(srGeometry.rotation*rad);
snr = sin(srGeometry.rotation*rad);
x   = xlc*csr - ylc*snr;
y   = xlc*snr + ylc*csr;

%  Convert from ENU to ECEF or UTM

if srGeometry.tf_latlon == 1
    
    %  ENU to ECEF
    
    alt     = 0;
    lat0rad = degtorad(srGeometry.latitude);
    lon0rad = degtorad(srGeometry.longitude);

    [x,y,z] = lv2ecef(x,y,zeros(size(x)),lat0rad,lon0rad,alt,srGeometry.ellipsoid);
    
    % ECEF to Geodetic.  h is not used but gives estimate of the effect of
    % Earth's curvature
    
    [mapy,mapx,h] = ecef2geodetic(x,y,z,srGeometry.ellipsoid);
    mapy          = radtodeg(mapy);
    mapx          = radtodeg(mapx);
    
elseif srGeometry.tf_latlon==0
    
    mapx = srGeometry.easting + x;
    mapy = srGeometry.northing  + y;
    
end

