function tlArrival = tlPick2tlArrival(srEvent, srStation, tlPickDir, stationIn, phaseIn, phaseOut, rlim, channelIn)
% Create a tlArrival structure from tlPick files
%
% Usage
%   tlArrival = tlPick2tlArrival(srEvent, srStation, tlPickDir, stationIn, phaseIn, phaseOut, rlim, channelIn)
% 
% Inputs
%   srEvent   - srEvent structure with derived x-y variable (use load_srEvent)
%   srStation - srStation structure with derived x-y variables (use load_srStation)
%   tlPickDir - directory with tlPick files
%   stationIn - cell array of station names (empty for all in srEvent)
%   phaseIn   - cell array of phases in tlPick files
%   phaseOut  - Corresponding phase names in tlArrival structure 
%               (empty implies phaseOut = phaseIn)
%   rlim      - Lower and upper range limits.  Options are
%               1 pair to use for all phases
%               A different pair for each phases in rows of a Nphase x 2 matrix
%               Empthy for no limits 
%   channelIn - Scalar input if tlArrival is to have picks for one specific channel
%               (empty implies no check of channelIn which means could get
%               more than one time for an arrival for channelIn specific picking)
%
% Outputs
%   tlArrival - tlArrival structure
%
% Uses inefficient concatenation of tlArrival vectors so may be slow for
% big data sets.
%
% EDITS (RTW 10/12/2011)
%       (1) Field "secs" to field "time" to be consistent with Stingray/TomoLab code.  
%       (2) "station" field now iterates & concatenates consistently with other fields.
%       (3) Phase verification now case-insensitive.
%
% Changed defalt behaviour to just ignore shots in the tlPick files that are not in srEvent - May 2013
%
% Updated to include picks only if the new flag '"use" is set to true in the tlPick file
%   Note: backward compatibility sets this flag to to true if not present in tlPick file

% Input arguments
if nargin<5
  error('tlPick2tlArrival requires at least 5 input arguments')
end
if nargin<6
  phaseOut = phaseIn;
end
if isempty(stationIn)
  stationIn = srStation.name;
end
if isempty(phaseOut)
  phaseOut = phaseIn;
end
if nargin<7
  rlim = [0 Inf];
end
if length(rlim(1,:))==2
  rlim = rlim(:)';
end
if nargin<8
  channelIn = [];
end

% Make string inputs cells
if ~iscell(stationIn)
  stationIn = {stationIn};
end
if ~iscell(phaseIn)
  phaseIn = {phaseIn};
end
if ~iscell(phaseOut)
  phaseOut = {phaseOut};
end

% Create empty tlArrival structure
tlArrival.eventtype = [];
tlArrival.eventid = [];
tlArrival.phase = [];
tlArrival.error = [];
tlArrival.station = [];
tlArrival.time = [];

% Only keep picks that are flagged use = true;
keep = [];

% Loop through phases and stations loading tlPick
for iPhase = 1:length(phaseIn)
  for iStation = 1:length(stationIn)
    if ~any(strcmp(srStation.name, stationIn{iStation}))
      error('tlPick2tlArrival - Invalid station (not in srStation)')
    end
    [station, channel, eventid, phase, time, unc, use] = ...
            read_tlPick(tlPickDir, stationIn{iStation}, phaseIn(iPhase), ...
             {'station', 'channelIn', 'eventid','phase','time','unc','use'});
    
    if ~isempty(station)
      % Verify tlPick station and phase       
      if ~all(strcmp(station, stationIn{iStation}))
        error('tlPick file with wrong station')
      end
      if ~all(strcmpi(phase, phaseIn{iPhase}))
        error('tlPick file with wrong phase')
      end

      % Select specified channelIn
      if ~isempty(channelIn)
        i = find(channel == channelIn);
        eventid = eventid(i);
        unc = unc(i);
        time = time(i);
        station = station(i);
        phase = phase(i);
        use = use(i);
      end

      % Verify tlPick events
      j = vector_indexmatch(srEvent.id, eventid);
      if any(j<0)
        error('Picks duplicated (possibly channelIn specific with channelIn = [])')
      end

%       if any(j==0)
%         error('tlPick includes an invalid event (not in srEvent)')
%       end
	  % Ignore events that are not in srEvent 
      if any(j==0)
        i=find(j);
        eventid = eventid(i);
        unc = unc(i);
        time = time(i);
        station = station(i);
        phase = phase(i);
        use = use(i);
        j = vector_indexmatch(srEvent.id, eventid);
      end
      
      % Select by range
      % Modified to allow for different ranges for different phases -
      % Arnoux, Jan 2017
      if ~isempty(rlim)
          if length(rlim(:)) > 2
              i = find(strcmp(srStation.name, stationIn{iStation}));
              range = sqrt((srEvent.x(j) - srStation.x(i)).^2 + (srEvent.y(j) - srStation.y(i)).^2);
              i = find(range>=rlim(1,iPhase) & range<=rlim(1,iPhase+length(rlim)/2));
              display(iStation); display(num2str(rlim(1,iPhase+length(rlim)/2)))
              eventid = eventid(i);
              unc = unc(i);
              time = time(i);
              station = station(i);
              phase = phase(i);
              use = use(i);
          else
              i = find(strcmp(srStation.name, stationIn{iStation}));
              range = sqrt((srEvent.x(j) - srStation.x(i)).^2 + (srEvent.y(j) - srStation.y(i)).^2);
              k = min(iPhase, size(rlim,1)); 
              i = find(range>=rlim(k,1) & range<=rlim(k,2));
              eventid = eventid(i);
              unc = unc(i);
              time = time(i);
              station = station(i);
              phase = phase(i);
              use = use(i);
          end
      end      

      % Add current tlPick data to output
      j = vector_indexmatch(srEvent.id, eventid);
      tlArrival.eventtype = [tlArrival.eventtype; srEvent.type(j)];
      tlArrival.eventid = [tlArrival.eventid; eventid];
      tlArrival.phase = [tlArrival.phase; repmat(phaseOut(iPhase),length(time),1)];
      tlArrival.error = [tlArrival.error; unc];
      tlArrival.station = [tlArrival.station; station];
      tlArrival.time = [tlArrival.time; time];
      keep = [keep; use];
    end
  end
end

tlArrival.num = length(tlArrival.time);
if any(~keep)
  tlArrival = subset_oneElementStructure(tlArrival, 'num', keep);
end
tlArrival.phaselist = unique(tlArrival.phase);
tlArrival.nphases = length(tlArrival.phaselist);
tlArrival.filename = '';