function srEvent=load_srEvent(theFile,srGeometry)

%LOAD_srEVENT -- Load Stingray structure srEvent
%                (Stingray utility)
%
%  srEvent = load_srEvent(theFile,srGeometry)
%
%  Load the srEvent structure; check required fields; fill derived fields.
%
%  INPUT:
%               theFile:    Filename
%            srGeometry:    Stingray structure
%
%  OUTPUT:
%               srEvent:    Stingray structure

%  Copyright 2010 Blue Tech Seismics, Inc.

load(theFile)

srEvent.filename = theFile;

%%  Required fields.  Not all are checked

if isempty(srEvent.id) || isempty(srEvent.type) 
    error('Required fields missing in srEvent')
end

%% Derived Fields

%  Check geographic vs. cartesian; get x y positions

if srGeometry.tf_latlon
    if isempty(srEvent.latitude) || isempty(srEvent.longitude)
        error('Required fields missing in srEvent')
    end
    [srEvent.x, srEvent.y] = map2xy(srEvent.longitude, srEvent.latitude, ...
        srGeometry);
elseif ~srGeometry.tf_latlon
    if isempty(srEvent.easting) || isempty(srEvent.northing)
        error('Required fields missing in srEvent')
    end

    [srEvent.x, srEvent.y] = map2xy(srEvent.easting, srEvent.northing, ...
        srGeometry);
elseif isempty(srGeometry.tf_latlon)
    error('ERROR:  srGeometry.tf_latlon has no value')
end

% Get model z-value of events

srEvent.z = nan(length(srEvent.id),1);

%  Airgun data

if any(srEvent.type == 1)
    I = srEvent.type==1;
    srEvent.z(I) = 0;
end

%  Borehole data

if any(srEvent.type == 2)
    I=find(srEvent.type==2);    
    srEvent.z(I)    = -srEvent.depth(I);
end

%  Regional data

if any(srEvent.type ==3)
    I=find(srEvent.type==3);    
    srEvent.z(I)    = -srEvent.depth(I);
end

%  Streamer data

if any(srEvent.type == 5)
    I = srEvent.type==5;
    srEvent.z(I) = 0;
end


% Number of events

srEvent.nevt = length(srEvent.id);

display(srEvent)

