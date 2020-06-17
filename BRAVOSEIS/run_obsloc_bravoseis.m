%% Script to run OBS location function for ETOMO experiment
% There is input structure p which contains stingray/tomolab structures and
% inputs that are documented below


%% Load Stringray and TomoLab inputs and put into p
% srControl is required for Stingray (only field "tf_latlon" is relevant)
srControl.tf_latlon = 1;
srControl.tf_anisotropy = 0;
srControl_tf_line_intergrate = 0;
srControl_arcfile = '';

% Load Geometry
% % srGeometry  = load_srGeometry('srGeometry');
% % %srGeometry  = mapproject(srGeometry,srControl);
% % 
% % % Load Stations and Events (Events)
% % srStation   = load_srStation('srStation',srGeometry);
% % srEvent     = load_srEvent('srEvent',srGeometry);
% % 
% % % Load Arrival Structure (water waves)
% % tlArrival   = load_tlArrival('tlArrival');

srElevation = fill_srElevation_fields(srElevation,srGeometry);

p.srControl = srControl;
p.srGeometry = srGeometry;
p.srStation = srStation;
p.srEvent = srEvent;
p.tlArrival = tlArrival;
p.srElevation = srElevation;

% What are we solving for
p.solve.xyStation = 1; 
p.solve.zStation = 1; 
p.solve.xyEvent = 0; 
p.solve.tEvent = 0; 
% Are we adjusting apriori depths to bathymetry?
p.solve.useBathymetry = 0; 
% Are we scaling the apriori data errors to match the misfits
p.solve.adjustDataError = 0;

% Apriori model parameter errors (assumed not to vary between stations/events)
p.Station.xErrorApr = 0.5;
p.Station.yErrorApr = 0.5;
p.Station.zErrorApr = 0.1;
p.Event.xErrorApr = 0.02;
p.Event.yErrorApr = 0.02;
p.Event.tErrorApr = 0.01;

% Travel time tables
p.nlookup = 1;
p.lookup(1).eventElevation = -0.015;
% % p.nlookup = 2;
addpath('/Users/earthnote/Desktop/Brainsfield/OBS_Event_Locations/OBS_Event_Relocation_Inversion/Ray_Tables')

p.lookup(1).file = 'TTtables_simple.mat';
% % p.lookup(1).file = 'ray_ttTable_AverageWaterVel1_9.mat'; 
% % p.lookup(1).eventElevation = -0.009;
% % p.lookup(2).file = 'ray_ttTable_AverageWaterVel1_15.mat' 
% % p.lookup(2).eventElevation = -0.015;

% Use LSQR when greater than this many model parameters
p.nlsqr = 11000;
% Maximum number of iterations
p.iterMax = 5;       
% Stop iterating when change in time and spatial model parameters is less than
p.dtMax = 1e-4;     
p.dxyzMax = 1e-4;
% Stop iterating when the change in rms and normalized rms is less than
p.drmsMax = 1e-4;
p.drmsNormMax = 1e-4;
% Do not use travel times for ranges greater than this
p.rangeMax = 12;
% Number of stations on residual versus shot index plots (0 for no plots)
p.nStaPlot = 100;
% Logical to add plot labels (slow)
p.plotLabel = 1;
% Controls vertical positon of plot labels
p.yLabelRms = 0.06;
p.yLabelRmsNorm = 10;

%% Run obsloc
s = obsloc(p);
s.reloc_description = 'Testing for Bransfield station relocation';

cd '/Users/earthnote/Desktop/Brainsfield/OBS_Event_Locations/OBS_Event_Relocation_Inversion/Obsloc_Structures';
save Obsloc_structure_stations_xy_events_xy s