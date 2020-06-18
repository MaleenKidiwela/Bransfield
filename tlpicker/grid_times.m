function [meanTime, nTime, stdTime] = grid_times(srEvent,srStation, ...
                           tMatrix,xGrid,yGrid,location,station,eventid)
% Grids a matrix of mean travel times by mid-point of shot location
% 
%
% Usage
%   [meanTime, nTime, stdTime] = grid_times(srEvent,srStation, ...
%                           tMatrix,xGrid,yGrid,location,station,eventid)
% 
% Inputs
%   srEvent   - Stingray srEvent structure
%   srStation - Stingray srStation structure
%   tMatrix   - Matrix of time with dimensions srEvent.nevt, srStation.nsta
%               NaN is used to indicate no data
%   xGrid     - X coordinates of center points of grid
%   ygrid     - Y coordinates of center points of grid
%   location  - Controls where times are assigned to
%               Either 'middle' for event-receiver midpoint (default)
%                 or 'event' for event location.
%   station   - Cell structure of stations to process (default is all)
%   eventid   - Vector of shot numbers to process (default is all).
%               Invalid shot numbers are ignored
%
% Outputs
%   meanTime  - Grid of mean times X in rows and Y in columns
%   nTime     - Number of times at each grid point
%   stdTime   - Grid of standard deviation of times.

%% Inputs
if nargin<6
  location = 'middle';
end
if nargin<7
  station = srStation.name;
end
if ischar(station)
  station = {station};
end
if nargin<8
  ievt = 1:srEvent.nevt;
else
  ievt = vector_indexmatch(srEvent.id, eventid);
  if any(~ievt);
    disp('grid_times eliminating invalid eventid values');
    eventid = eventid(~~ievt);
    ievt = ievt(~~ievt);
  end
end

%% Empty outputs
meanTime = zeros(length(xGrid),length(yGrid));
nTime = zeros(length(xGrid),length(yGrid));
stdTime = zeros(length(xGrid),length(yGrid));

%% Loop through stations
for iloop = 1:length(station)
  ista = find(strcmp(srStation.name, station{iloop}));
  if isempty(ista)
    error('Invalid station name');
  end
  
  % Only process non NaNs
  index = find(~isnan(tMatrix(ievt,ista)));
  
  % X, Y positions for times
  if strcmpi(location,'middle') || strcmpi(location,'midpoint')
    x = (srEvent.x(ievt(index)) + srStation.x(ista))/2;
    y = (srEvent.y(ievt(index)) + srStation.y(ista))/2;
  elseif strcmpi(location,'event') || strcmpi(location,'shot');
    x = srEvent.x(ievt(index));
    y = srEvent.y(ievt(index));
  end
  
  % Index of closest grid point
  i = dsearchn(xGrid(:),x);
  j = dsearchn(yGrid(:),y);
  
  % Assign times
  for k = 1:length(i);
    meanTime(i(k),j(k)) = meanTime(i(k),j(k)) + tMatrix(ievt(index(k)),ista);
    nTime(i(k),j(k)) = nTime(i(k),j(k)) + 1;
    stdTime(i(k),j(k)) = stdTime(i(k),j(k)) + tMatrix(ievt(index(k)),ista)^2;
  end
end

meanTime = meanTime ./ nTime;
meanTime(~nTime) = NaN;
stdTime = sqrt(stdTime./nTime - meanTime.^2);
stdTime(~nTime) = NaN;
    
    
    