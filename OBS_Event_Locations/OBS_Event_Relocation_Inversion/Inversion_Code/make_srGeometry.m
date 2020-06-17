% make_srGeometry
% specific for the ETOMO experiment based on the locations of the outer
% OBS lines.  
% Geographic Coordinate System is spherical
% 			
% srGeometry.longitude (decimal degrees)
%       Chosen longitude for center of experiment.  Becomes
%       origin for x-axis in cartesian coordinate system. Value does not
%       affect short distance conversion factors.  Choose one that is
%       convenient.
% 
% srGeometry.latitude (decimal degrees )
%       Chosen latitude for center of experiment.  Becomes
%       origin for y-axis in cartesian coordinate system. Value is
%       important to calculation of short distance conversion factors.
%       Choose one that is near to geographic center of experiment; this
%       will reduce error of mapping between geographic and cartesian
%       coordinates.
% 
% srGeometry.rotation (degrees )
%       Rotation of cartesian coordinate system with respect to
%       geographic coordinate system.  When the rotation is zero, the
%       x-axis points east and the y-axis points north. A positive angle
%       rotates the x-y axis in the CCW direction.   This is useful for
%       limiting the dimensions of the cartesian space that is used for ray
%       tracing.  See Coordinate System

%% from cruiseplan11.m  
% Coordinate paramters - Outer (undershoot) and inner can be different
% at MEF origin_outer.lat = 47.9489;
%        origin_outer.lon = -129.0987;
% at Salty Dawg on MCS 3
origin_outer.lat = 47.978;
origin_outer.lon = -129.075;
origin_outer.rot = -17.5;
origin_inner = origin_outer;
origin_outer.rot = -5.5;
%
% OBS locations and type
% % Specify x & y(km) with outer lines come first
nobs = 64;
nobs_outer = 28;
 yobs = [[-50:10:40]+5 [-30 -10 8 30] -30:20:10 28 [-50 -42 -30:10:20 28 44]+5 -30:7.5:22.5 [-10 -5 0.5 5 10] -26.25 -18.75 -12.5:5:12.5 18.75 26.26 -10:5:10 -22.5:7.5:22.5 29 ];
 xobs = [-23*ones(1,10) -16*ones(1,4) 16*ones(1,4) 23*ones(1,8) 25 18 [-9 -9 -10 -8 -9 -7.5 -9 -9] -4*ones(1,5) 0*ones(1,10) 4*ones(1,5) 8*ones(1,7) 6 ];

figure(1),clf
plot(xobs,yobs,'*')
grid on

nline = 35;%total number of lines without the MCS line
nline_outer = 6;%outer lines
nshot_outer_st = 21;
nshot_outer_end = 259;

nline_middle = 19;%middle lines
nshot_middle_st = 21;
nshot_middle_end = 156;

nline_inner = 10;%inner lines
nshot_inner_st = 66;
nshot_inner_end = 111;

%Specifying undershoot shot lines; converted in 'lonline' and 'latline'  
yline = [-50*ones(1,6) -30*ones(1,19) -10*ones(1,10);
         50*ones(1,5) 53 30*ones(1,19) 10*ones(1,10)];
xline = [-30 -23 -16 16 23 30 -9:1:9 -4.5:1:4.5;
         -30 -23 -16 16 23 30 -9:1:9 -4.5:1:4.5];

%Specifying shot locations in the outer area;
ybase_outer = [-49.96:0.420:50];
xbase_outer = [-30 -23 -16 16 23 30];
ybase_outer2 = ybase_outer.'

youter = repmat(ybase_outer2, 1, nline_outer)
xouter = repmat(xbase_outer, nshot_outer_end - nshot_outer_st + 1,1)


%Specifying shot locations in the middle area;
ybase_middle = [-30.75:0.450:30];
xbase_middle = [-9:1:9];
ybase_middle2 = ybase_middle.'

ymiddle = repmat(ybase_middle2, 1, nline_middle)
xmiddle = repmat(xbase_middle, nshot_middle_end - nshot_middle_st + 1,1)


%Specifying shot locations in the inner area;
ybase_inner = [-10.25:0.45:10];
xbase_inner = [-4.5:1:4.5];
ybase_inner2 = ybase_inner.'

yinner = repmat(ybase_inner2, 1, nline_inner)
xinner = repmat(xbase_inner, nshot_inner_end - nshot_inner_st + 1,1)


% % MCS line 3 specified in polyfit/polyval format
% [x,y] = map2xy(lon_mcs3,lat_mcs3,inner);
% p = polyfit(x,y,1);
% x = [-0.5 0.5] * length_mcs3;
% y = polyval(p,x);
% xline = [xline x'];
% yline = [yline y'];

figure(1), hold on, plot(xline,yline)
xlim([-35 35])

% Coordinate transformation 
%  Separate for outer (undershoot) and inner portions but ONLY VARY ROTATION
[xltkm_outer,xlnkm_outer] = setorg(origin_outer.lat);
outer = [origin_outer.lon origin_outer.lat xlnkm_outer xltkm_outer origin_outer.rot];
[xltkm_inner,xlnkm_inner] = setorg(origin_inner.lat);
inner = [origin_inner.lon origin_inner.lat xlnkm_inner xltkm_inner origin_inner.rot];
xltkm = xltkm_outer;
xlnkm = xlnkm_outer;

% Convert to latitude/longitude
[lonobs(1:nobs_outer),latobs(1:nobs_outer)] = xy2map(xobs(1:nobs_outer),yobs(1:nobs_outer),outer);
[lonobs(nobs_outer+1:nobs),latobs(nobs_outer+1:nobs)] = xy2map(xobs(nobs_outer+1:nobs),yobs(nobs_outer+1:nobs),inner);
[lonline(:,1:nline_outer),latline(:,1:nline_outer)] = xy2map(xline(:,1:nline_outer),yline(:,1:nline_outer),outer);
[lonline(:,nline_outer+1:nline),latline(:,nline_outer+1:nline)] = xy2map(xline(:,nline_outer+1:nline),yline(:,nline_outer+1:nline),inner);
[lonsh_outer(:,1:nline_outer),latsh_outer(:,1:nline_outer)] = xy2map(xouter(:,1:nline_outer),youter(:,1:nline_outer),outer);
[lonsh_middle(:,1:nline_middle),latsh_middle(:,1:nline_middle)] = xy2map(xmiddle(:,1:nline_middle),ymiddle(:,1:nline_middle),inner);
[lonsh_inner(:,1:nline_inner),latsh_inner(:,1:nline_inner)] = xy2map(xinner(:,1:nline_inner),yinner(:,1:nline_inner),inner);

figure(2),clf
plot(lonobs,latobs,'o','markerfacecolor','r')
hold on
plot(lonsh_outer,latsh_outer,'r')
plot(lonsh_middle,latsh_middle,'g')
plot(lonsh_inner,latsh_inner,'b')

%% Check compared to the real lines
load srEvent
load srStation

figure(2)
plot(srStation.longitude,srStation.latitude,'o','markerfacecolor','k')
plot(srEvent.longitude,srEvent.latitude,'.k')

%%%%% STILL TO DO:  Use map2xy to plot the true locations on the x-y grid. %%%%%

%% make and save srGeometry
srGeometry.longitude = origin_outer.lon;
srGeometry.latitude = origin_outer.lat ;
srGeometry.rotation = origin_outer.rot;

save srGeometry srGeometry

