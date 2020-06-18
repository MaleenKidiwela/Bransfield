function [hFigMap, haxes] = plot_shotDataMap(srStation, srEvent, data, station, varargin)
% Plots color coded data at shot source locations or source/receiver mid points
%
% Usage
%   [hFigMap, haxes] = plot_shotDataMap(srStation, srEvent, data, station, prop1, val1, prop2, val2, ...)
%
% Inputs
%   srStation - Stingray Station structure (see Stingray manual)
%   srEvent   - Stingray Event structure (see Stingray manual)
%   data      - Vector of values for events (length is srStation.nevt)
%               NaN if no data for a particular event
%   station   - Station name (string)
%   prop, val - string-values pairs setting parameters to non-default values
%               PROP may be the following case insensitive string
%               'lim'  - VAL is the [minimum maximum] of color scale
%                        Default is the minimum / maximum of DATA
%               'outOfBounds' - VAL is a string that indicates how out of bounds
%                          data are treated
%                        Default is 'error' to exit on an error
%                        'fix' - fixes out of bounds values to limits
%                        'warn'   - warns of out of bounds values and fixes
%               'cmap' - VAL is the choice of color map
%                        Default is 'jet'
%               'midPoint' - VAL is logical to plot data at
%                            midPoint rather than event location.
%                            Default is false
%               'markerSize' - VAL is the size of data circles (Default is 8)
%               'noDataMarkerSize' - VAL is the size of no data points (Default is 5)
%               'noDataMarkerColor' - VAL is the color of no data points (Default is 5)
%               'stationMarkerSize' - VAL is the size of station squares (Default is 5)
%               'stationMarkerColor' - VAL is the color of station squares (Default is 5)
%               'usedStationMarkerSize' - VAL is the size of station squares (Default is 5)
%               'usedStationMarkerColor' - VAL is the color of station squares (Default is 5)
%               'plotTitle' - VAL is plot title string (Default is '');)
%               'cbarTitle' - VAL is colorbar title string (Default is '');
%
% Outputs
%   hFigMap - Handle of figure
%

%% Input argument processing
if nargin < 4
  error('plot_shotDataMap requires at least 4 arguments');
end
ista = find(strcmpi(srStation.name, station));
lim = [min(data) max(data)];
outOfBounds = 'error';
cmap = 'jet';
midPoint = false;
markerSize = 8;
noDataMarkerSize = 5;
noDataMarkerColor = 'k';
stationMarkerSize = 5;
stationMarkerColor = 'k';
usedStationMarkerSize = 5;
usedStationMarkerColor = 'k';
plotTitle = '';
cbarTitle = '';

for i = 1:2:nargin-5
  if strcmpi(varargin(i),'lim')
    lim = varargin{i+1};
  elseif strcmpi(varargin(i),'cmap')
    cmap = varargin{i+1};
  elseif strcmpi(varargin(i),'outofbounds')
    outOfBounds = varargin{i+1};
  elseif strcmpi(varargin(i),'midPoint')
    midPoint = varargin{i+1};
  elseif strcmpi(varargin(i),'markerSize')
    markerSize = varargin{i+1};
  elseif strcmpi(varargin(i),'noDataMarkerSize')
    noDataMarkerSize = varargin{i+1};
  elseif strcmpi(varargin(i),'noDataMarkerColor')
    noDataMarkerColor = varargin{i+1};
  elseif strcmpi(varargin(i),'stationMarkerSize')
    stationMarkerSize = varargin{i+1};
  elseif strcmpi(varargin(i),'stationMarkerColor')
    stationMarkerColor = varargin{i+1};
  elseif strcmpi(varargin(i),'usedStationMarkerSize')
    stationMarkerSize = varargin{i+1};
  elseif strcmpi(varargin(i),'usedStationMarkerColor')
    stationMarkerColor = varargin{i+1};
  elseif strcmpi(varargin(i),'plotTitle')
    plotTitle = varargin{i+1};
  elseif strcmpi(varargin(i),'cbarTitle')
    cbarTitle = varargin{i+1};
  end  
end

% Color values for data and out of bounds testing
index = find(data<lim(1));
if ~isempty(index)
  if strcmpi(outOfBounds,'error') || strcmpi(outOfBounds,'warn')
    fprintf('Events below lower bound\n')
    for i=1:length(index)
      fprintf('Event %6.0f with value %10.4f\n',srEvent.id(index(i)),data(index(i)))
    end
  end
  if strcmpi(outOfBounds,'error')
    error('Data values below lower limits')
  end
  if strcmpi(outOfBounds,'fix') || strcmpi(outOfBounds,'warn')
    data(index) = lim(1);
  end
end
index = find(data>lim(2));
if ~isempty(index)
  if strcmpi(outOfBounds,'error') || strcmpi(outOfBounds,'warn')
    fprintf('Events above upper bound\n')
    for i=1:length(index)
      fprintf('Event %6.0f with value %10.4f\n',srEvent.id(index(i)),data(index(i)))
    end
  end
  if strcmpi(outOfBounds,'error')
    error('Data values above lower limits')
  end
  if strcmpi(outOfBounds,'fix') || strcmpi(outOfBounds,'warn')
    data(index) = lim(2);
  end
end
col = LinColor(data,lim,cmap);

% Plotting
if ~midPoint
  x = srEvent.x;
  y = srEvent.y;
else
  x = (srEvent.x + srStation.x(ista))/2;
  y = (srEvent.y + srStation.y(ista))/2;
end
plot(x, y, '.k', 'markerSize', noDataMarkerSize, ...
     'markeredgecolor',noDataMarkerColor,'markerfacecolor',noDataMarkerColor);
haxes = gca;

hold on
plot(srStation.x, srStation.y, 'sk', 'markerSize', stationMarkerSize, ...
     'markeredgecolor',stationMarkerColor,'markerfacecolor','none');
for i = 1:length(data)
  if ~isnan(data(i))
    plot(x(i), y(i), 'o','markerfacecolor', col(i,:), 'markeredgecolor', col(i,:), ...
         'markerSize',markerSize);
  end
end
plot(srStation.x(ista), srStation.y(ista), 'pk', 'markerSize', usedStationMarkerSize, ...
     'markeredgecolor',usedStationMarkerColor,'markerfacecolor',usedStationMarkerColor);

% Labeling and colorbar
axis('equal');
xlabel('X, km');
ylabel('Y, km');
title(plotTitle);
caxis(lim);
hcb = colorbar;
h = gca;
set(gcf,'currentaxes',hcb);
title(cbarTitle);
set(gcf,'currentaxes',h)

hFigMap = gcf;

