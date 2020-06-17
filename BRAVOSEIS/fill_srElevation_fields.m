function [srElevation] = fill_srElevation_fields(srElevation,srGeometry)

%FILL_srELEVATION -- Fills derived fields of the srElevation structure
%                    (Stingray utility)
%
%  Required fields:
%
%         srElevation.header(1:8)
%         srElevation.data
%
%  Derived fields:
%
%  If srGeometry.tf_latlon=true;
%
%         srElevation.longitude   (vector)
%         srElevation.latitude    (vector)
%         srElevation.LON         (mesh)
%         srElevation.LAT         (mesh)
%         srElevation.X           (mesh)
%         srElevation.Y           (mesh)
%
%  If srGeometry.tf_latlon=false;
%
%         srElevation.easting     (vector)
%         srElevation.northing    (vector)
%         srElevation.EASTING     (mesh)
%         srElevation.NORTHING    (mesh)
%         srElevation.X           (mesh)
%         srElevation.Y           (mesh)


if srGeometry.tf_latlon == 1
    
    % x-direction.  Row vector
    
    srElevation.longitude   = linspace(srElevation.header(1), ...
        srElevation.header(1)+(srElevation.header(7)-1)*srElevation.header(5),...
        srElevation.header(7));
    
    % y-direction.  Transpose to be a column vector
    
    srElevation.latitude    = linspace(srElevation.header(3), ...
        srElevation.header(3)+(srElevation.header(8)-1)*srElevation.header(6),...
        srElevation.header(8))';
    
    % mesh, switch xy calls on meshgrid since data is nx by ny
    
    [srElevation.LAT srElevation.LON] = meshgrid(srElevation.latitude,srElevation.longitude);
    [srElevation.X   srElevation.Y]   = map2xy(srElevation.LON,srElevation.LAT,srGeometry);
    
    
elseif srGeometry.tf_latlon == 0
    
    srElevation.easting   = linspace(srElevation.header(1), ...
        srElevation.header(1)+(srElevation.header(7)-1)*srElevation.header(5),...
        srElevation.header(7));
    
    % y-direction.  Transpose to be a column vector
    
    srElevation.northing    = linspace(srElevation.header(3), ...
        srElevation.header(3)+(srElevation.header(8)-1)*srElevation.header(6),...
        srElevation.header(8))';
    
    % mesh, switch xy calls on meshgrid since data is nx by ny
    
%     display('Check transpose of srElevationn')
%     keyboard
    [srElevation.NORTHING srElevation.EASTING] = meshgrid(srElevation.northing,srElevation.easting);
    [srElevation.X   srElevation.Y]   = map2xy(srElevation.EASTING,srElevation.NORTHING,srGeometry);
    
end


