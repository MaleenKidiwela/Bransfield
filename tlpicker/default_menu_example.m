% Sets default menu defaultMenueters

% Values for Orienting Channels
clear defaultMenu

% Experiment
defaultMenu.experiment.srGeometry = '/users/wilcock/dropbox/research/ETOMO/etomo/srGeometry_ETOMO_Nov2011.mat';
defaultMenu.experiment.srStation = '/users/wilcock/dropbox/research/ETOMO/etomo/srStation_ETOMO_relocated.mat';
defaultMenu.experiment.srEvent = '/users/wilcock/dropbox/research/ETOMO/etomo/srEvent_ETOMO_relocated_stacked.mat';
defaultMenu.experiment.segyCatalog = '/Users/wilcock/Dropbox/Software/tlpicker/catalog/segyCatalog_etomo_Nov2011.mat';

% Trace Selection
defaultMenu.segy.station = 'BRA13';
defaultMenu.segy.eventid = '23000:23999';
defaultMenu.segy.channel = '3';
defaultMenu.segy.excludeBad = 1;
defaultMenu.segy.tlim0 = [ ]; 
defaultMenu.segy.tlim1 = [ ];
defaultMenu.segy.redVel = [ ];
defaultMenu.segy.select(1).name = [ ];
defaultMenu.segy.select(1).lim0 = [ ];
defaultMenu.segy.select(1).lim1 = [ ];
defaultMenu.segy.select(2).name = [ ];
defaultMenu.segy.select(2).lim0 = [ ];
defaultMenu.segy.select(2).lim1 = [ ];
defaultMenu.segy.select(3).name = [ ];
defaultMenu.segy.select(3).lim0 = [ ];
defaultMenu.segy.select(3).lim1 = [ ];
defaultMenu.segy.select(4).name = [ ];
defaultMenu.segy.select(4).lim0 = [ ];
defaultMenu.segy.select(4).lim1 = [ ];
defaultMenu.segy.pickSelect.DFS = '';
defaultMenu.segy.pickSelect.phase = '';
defaultMenu.segy.pickSelect.chanIter = 0;

% Filtering
defaultMenu.filter.lim0 = 5;
defaultMenu.filter.lim1 = 40;
defaultMenu.filter.order = 4;
defaultMenu.filter.muteTime = 0.1;
defaultMenu.filter.zeroPhase  = false;
defaultMenu.filter.onLoad = true;  

% Statics and Muting
defaultMenu.static.DFS = '';
defaultMenu.static.phase = '';
defaultMenu.static.chanIter = 0;
defaultMenu.mute.DFS = [];
defaultMenu.mute.phase = [];
defaultMenu.mute.chanIter = 0;

% Sorting 
defaultMenu.sort(1).name = '';
defaultMenu.sort(2).name = '';
defaultMenu.sort(3).name = '';
defaultMenu.sort(4).name = '';

% Plotting
defaultMenu.plot.xOption = 'xevt';
defaultMenu.plot.xlim0 = 0;
defaultMenu.plot.xlim1 = 0;
defaultMenu.plot.xScaleOption = 'x';
defaultMenu.plot.xScale = 1e-5;
defaultMenu.plot.clip = 0.1;
defaultMenu.plot.wiggleOption = 1;
defaultMenu.plot.tlim0 = 0;
defaultMenu.plot.tlim1 = 20;
defaultMenu.plot.redVel = 0;
defaultMenu.plot.defaultTitle = true;
defaultMenu.plot.demedian = false;
defaultMenu.plot.title= '';
defaultMenu.plot.labelIncrement = 3;
defaultMenu.plot.labelOption = 'range';

% Inactive arrivals (picks)
defaultMenu.inactive.DFS = '/users/wilcock/dropbox/research/ETOMO/etomo/Picks_syn';
defaultMenu.inactive.phase = '';
defaultMenu.inactive.chanIter = [];
defaultMenu.inactive.lineType = 'b g';

% Postscript file
defaultMenu.ps.fileName = '';
defaultMenu.ps.append = false;

% Active arrivals and picking
defaultMenu.active.directory = '/users/wilcock/dropbox/research/ETOMO/etomo/Picks';
defaultMenu.active.phase = 'Pw';
defaultMenu.active.channelSpecific = 0;
defaultMenu.active.unc = 0.1;
[~,user] = unix('whoami');
defaultMenu.active.user = user(1:end-1);
defaultMenu.active.comment = '';
defaultMenu.active.orderEventID = true;

% Views
defaultMenu.view.fileName = 'myview';

% Experiment Geometry Plot
defaultMenu.geometryPlot.loaded = true;
defaultMenu.geometryPlot.inactive = false;
defaultMenu.geometryPlot.active = false;

% New Stuff - Not yet in the menu
defaultMenu.orientHoriz.srStationRotation = ' ';
defaultMenu.orientHoriz.tlim0 = 0;
defaultMenu.orientHoriz.tlim1 = 0.05;
defaultMenu.orientHoriz.minEvt = 10;
defaultMenu.orientHoriz.minScale = 10^-1;
defaultMenu.orientHoriz.maxScale = 10^1;
defaultMenu.orientHoriz.nScale = 201;
defaultMenu.orientHoriz.rangeCrit = 2;
defaultMenu.orientHoriz.plotEachEvt = false;
defaultMenu.orientHoriz.noVertical = false;



%% Basic values
% % clear defaultMenu
% % 
% % % Experiment
% % defaultMenu.experiment.srGeometry = '/users/matlab/etomo/srGeometry_ETOMO.mat';
% % defaultMenu.experiment.srStation = '/users/matlab/etomo/srStation_ETOMO_relocated.mat';
% % defaultMenu.experiment.srEvent = '/users/matlab/etomo/srEvent_ETOMO_relocated.mat';
% % defaultMenu.experiment.segyCatalog = '/users/matlab/tlpicker/segyCatalog_etomo';
% % 
% % % Trace Selection
% % defaultMenu.segy.station = '21';
% % defaultMenu.segy.eventid = '1000:3999';
% % defaultMenu.segy.channel = '1';
% % defaultMenu.segy.excludeBad = 1;
% % defaultMenu.segy.tlim0 = 0; 
% % defaultMenu.segy.tlim1 = 5;
% % defaultMenu.segy.redVel = 6;
% % defaultMenu.segy.select(1).name = 'range';
% % defaultMenu.segy.select(1).lim0 = 0;
% % defaultMenu.segy.select(1).lim1 = 30;
% % defaultMenu.segy.select(2).name = [ ];
% % defaultMenu.segy.select(2).lim0 = [ ];
% % defaultMenu.segy.select(2).lim1 = [ ];
% % defaultMenu.segy.select(3).name = [ ];
% % defaultMenu.segy.select(3).lim0 = [ ];
% % defaultMenu.segy.select(3).lim1 = [ ];
% % defaultMenu.segy.select(4).name = [ ];
% % defaultMenu.segy.select(4).lim0 = [ ];
% % defaultMenu.segy.select(4).lim1 = [ ];
% % defaultMenu.segy.pickSelect.DFS = [];
% % defaultMenu.segy.pickSelect.phase = [];
% % defaultMenu.segy.pickSelect.chanIter = [];
% % 
% % % Filtering
% % defaultMenu.filter.lim0 = 5;
% % defaultMenu.filter.lim1 = 40;
% % defaultMenu.filter.order = 4;
% % defaultMenu.filter.muteTime = 0.2;
% % defaultMenu.filter.zeroPhase  = false;
% % defaultMenu.filter.onLoad = true;  
% % 
% % % Statics and Muting
% % defaultMenu.static.DFS = [];
% % defaultMenu.static.phase = [];
% % defaultMenu.static.chanIter = 0;
% % defaultMenu.mute.DFS = [];
% % defaultMenu.mute.phase = [];
% % defaultMenu.mute.chanIter = 0;
% % 
% % % Sorting 
% % defaultMenu.sort(1).name = 'range';
% % defaultMenu.sort(2).name = '';
% % defaultMenu.sort(3).name = '';
% % defaultMenu.sort(4).name = '';
% % 
% % % Plotting
% % defaultMenu.plot.xOption = 'sort';
% % defaultMenu.plot.xlim0 = 0;
% % defaultMenu.plot.xlim1 = 0;
% % defaultMenu.plot.xScaleOption = 'f';
% % defaultMenu.plot.xScale = 0.001;
% % defaultMenu.plot.clip = 2;
% % defaultMenu.plot.wiggleOption = 1;
% % defaultMenu.plot.tlim0 = 1;
% % defaultMenu.plot.tlim1 = 4;
% % defaultMenu.plot.redVel = 6;
% % defaultMenu.plot.defaultTitle = true;
% % defaultMenu.plot.demedian = false;
% % defaultMenu.plot.title= '';
% % defaultMenu.plot.labelIncrement = 10;
% % defaultMenu.plot.labelOption = 'eventid xmid';
% % 
% % % Inactive arrivals (picks)
% % defaultMenu.inactive.DFS = 'testdata/upicker';
% % defaultMenu.inactive.phase = 'pg';
% % defaultMenu.inactive.chanIter = [];
% % defaultMenu.inactive.lineType = 'b';
% % 
% % % Postscript file
% % defaultMenu.ps.fileName = '';
% % defaultMenu.ps.append = false;
% % 
% % % Active arrivals and picking
% % defaultMenu.active.directory = 'testdata/upicker';
% % defaultMenu.active.phase = 'dpg';
% % defaultMenu.active.channelSpecific = 0;
% % defaultMenu.active.unc = 0.1;
% % [~,user] = unix('whoami');
% % defaultMenu.active.user = user(1:end-1);
% % defaultMenu.active.comment = '';
% % defaultMenu.active.orderEventID = true;
% % 
% % % Views
% % defaultMenu.view.fileName = 'myview';
% % 
% % % Experiment Geometry Plot
% % defaultMenu.geometryPlot.loaded = true;
% % defaultMenu.geometryPlot.inactive = true;
% % defaultMenu.geometryPlot.active = true;




